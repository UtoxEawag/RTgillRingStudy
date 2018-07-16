##########################################################
# DRC curves with profile likelihood confidence intervals
#   
# Authors:                                              #
# Diana Coman Schmid                                    #
# Anze Zupanic
#
#  IN THE CODE YOU WILL FIND AT THE END OF LINES, SEVERAL SYMBOLS
# *** - indicates that the user should change the code before running it
# !!! - indicates possible troubleshooting lines - places where the code can fail, and reasons why
##########################################################


# INSTRUCTIONS 
# PUT THIS SCRIPT AND ALL THE DATA THAT YOU WANT ANALYZED INTO A FOLDER
# SEE DRC_BATCH_INPUTFILE.XLSX AND FORMAT YOUR DATA ACCORDINGLY
# THIS IS SCRIPT FOR PROCESSING MANY SUCH FILES SEQUENTIALLY

# clear the workspace
rm(list = ls())

# load libraries
if ( ! require(readr) )     { install.packages("readr");     library(readr) }
if ( ! require(readxl) )     { install.packages("readxl");     library(readxl) }
if ( ! require(drc) )     { install.packages("drc");     library(drc) }
if ( ! require(nlmrt) )     { install.packages("nlmrt");     library(nlmrt) }
if ( ! require(nls2) )     { install.packages("nls2");     library(nls2) }
if ( ! require(ggplot2) )     { install.packages("ggplot2");     library(ggplot2) }
if ( ! require(reshape) )     { install.packages("reshape");     library(reshape) }
if ( ! require(reshape2) )     { install.packages("reshape2");     library(reshape2) }
if ( ! require(XLConnect) )     { install.packages("XLConnect");     library(XLConnect) }
if ( ! require(lattice) )     { install.packages("lattice");     library(lattice) }
if ( ! require(sfsmisc) )     { install.packages("sfsmisc");     library(sfsmisc) }


# ******************** USER INPUT *****************************
# make working folders before running the script !!!
# define the working directories
# when copy file path from Windows Explorer please replace \ with //
dat.path <- "\\\\eawag\userdata//zupanian//My Documents//anze//scripts//R scripts//dose_response_batch"
res.path <- "\\\\eawag\userdata//zupanian//My Documents//anze//scripts//R scripts//dose_response_batch//output" 
# your directory where the output results will be stored


## read all files
# collect the file names available in the working directory
# !!!!!!!!  Close all excel files before running this code 
list.files(dat.path)
file.list = list.files(dat.path, pattern = "*.xlsx")
names(file.list) <- file.list
# as a results of this you get a list with all file names

## format all files
# read in the excel workbooks (files) containing sheets for each chemical/lab
dat.in <- list()
# for each file in the folder
#   load the excel file
#   read each sheet
#   write data into dat.in
for (f in names(file.list)){            
  wb <- loadWorkbook(file.path(dat.path,f))
  ws = readWorksheet(wb, sheet = getSheets(wb))  # get all sheets at once
  # in case there is a single sheet
  if (is.data.frame(ws)) {
    nameSheet = getSheets(wb)
    ws2 = list(ws)
    names(ws2) = nameSheet
    ws = ws2
  }
  dat.in[[f]] <- ws
}
# as a results of this you get another list (dat.in)
# in this list there one list for each file (with the same number of dataframes as there are sheets in the file)

## remove spaces from excel workbook (files) and worksheet names  
names(dat.in) <- gsub(" ","_",names(dat.in))
for (n in names(dat.in)){
  names(dat.in[[n]]) <- gsub(" ","_",names(dat.in[[n]]))
}

# make an output dataset
dat.out <- list()

# for each file repeat the following analysis
for (b in names(dat.in)){
  # for each file make own output subdataset
  dat.out.t <- list()
  # open a new excel file (this is where the results of the analysis will go)
  res.wb <- loadWorkbook(file.path(res.path,paste("Res_",b,sep="")),create=TRUE)
  # on screen output tells which file is processed
  cat(paste("Processing Workbook >>", b, "\n", sep=''))
  # for each sheet in the file
  for (s in names(dat.in[[b]])){
    # on screen output tells which sheet is processed
    cat(paste("Processing Workbook >>", b,"  Worksheet >>",s, "\n", sep=''))
    # this is now a dataframe with the data of the sheet
    df <- dat.in[[b]][[s]][,-2]# omit the "stdev.conc" column
    df.long <- na.omit(melt(df,id="conc"))# switch from wide to long format and omit "NA" entries
    # omit any line if the concentration is 0 (control)
    if (any(df.long$conc == 0)) { # if there is a control (0) concentration, then it is removed
      index0 = which(df.long$conc == 0)
      df.long = df.long[-index0,]
    }
    ### !!!!!!!!!!!!!! this only works if the header names in the files are correct, so check the example_input_file carefully!!!
    rep <- max(  substring(as.character(df.long$variable), 2, 2)) #get the number of bioReps are available    
    df.long$bioRep <- as.factor(paste("bioRep",substring(as.character(df.long$variable), 2, 2),sep="")) #add column with biological replication info
    colnames(df.long) <-c("conc","var","pViab","bioRep") #change column names in the dataframe
    # log the concentration values
    df.long$conc = log10(df.long$conc)
    # working with biological replicates
    rep.n <- unique(df.long$bioRep)
    names(rep.n) <- unique(df.long$bioRep)
    ## for each biological replicate separately
    for (br in names(rep.n)){
      # on screen output tells which biological replicate is beeing processes
      cat(paste("Calculating 95% Profile Likelihood CI for Workbook >>", b,"  Worksheet >>",s,"Replicate",br, "\n", sep=''))
      # EC50 can only beetween highest and lowest exposure concentrations
      conc_opt <- seq(min(df.long$conc),max(df.long$conc),length.out=1000)  
      # initialiazing goodness of fit and slope values (see below the nonlinear drc equation)
      goodnessfit <- rep(NA, length(conc_opt))
      slope_value <- rep(NA, length(conc_opt))
      data_for_fit = df.long[which(df.long$bioRep == br),]
      fitted_models <- list()
      # model to be fitted
      drc_formula = pViab ~ 100/(1+10^((fixEC50-conc)*slope))
      # fit the curve to find the best fitting EC50. this is done by ftting the slope for different fixed values of EC50 (many optimizations)
      for(i in 1:length(conc_opt)) {
        fixEC50 <- conc_opt[i]
        tryCatch({
          # a combination of the nlxb and nls2 function is used for the fitting to guarantee a fit is found (99.8% success in finding a fit for different datasets so far)
          fit1 = nlxb(drc_formula,
                      start = c(slope = -4),     
                      trace = FALSE,
                      data = data_for_fit,
                      weights = c(1/pViab))
          mod1.pl <- nls2(drc_formula, data = data_for_fit, start = fit1$coefficients, weights = c(1/pViab),
                          algorithm = "brute-force")
		  fitted_models[[i]] = mod1.pl
        },error=function(e){cat("Error :", conditionMessage(e),"\n")} )
        goodnessfit[i] <- sqrt(sum(residuals(mod1.pl)^2))
        coefficient = coef(mod1.pl)
        slope_value[i] = coefficient[1]
      }
      index1 <- which(goodnessfit<(min(goodnessfit)+qchisq(p=0.95, df=3))) # calculate the profile likelihood CIs (95%)
      index2 <- which(goodnessfit==(min(goodnessfit)))# calculate the best fit
      # 95% CI lower and uppler value
      best_EC50_lower = conc_opt[min(index1)]
      best_EC50_upper = conc_opt[max(index1)]
      # best fitting EC50
      best_EC50 = conc_opt[median(index2)]
      # best fitting slope
      best_slope = slope_value[median(index2)]
      # quality of fit
      min(goodnessfit)
	    fixEC50 = best_EC50
	    # best fitting model
      curvefit= fitted_models[[median(index2)]]
      
      # output summary of the best fitting model 
      dat.out.t[[s]] <- summary(curvefit)
      
      # some outputs for quality control etc, optional
      
      # save results into excel file (EC50s, EC50 CIs, Slope)
      cat(paste("Saving model estimates to Workbook >>", b,"  Worksheet >>",s, "\n", sep=''))
      pl_ci <- as.data.frame(rbind(cbind(paste("EC50_PL95%CI_Lower Limit:",br,sep=""),as.character(best_EC50_lower)),
                                   cbind(paste("EC50_PL95%CI_Upper Limit:",br,sep=""),as.character(best_EC50_upper))))
      pl_ci <- as.data.frame(rbind(pl_ci,cbind(paste("EC50_PL:",br,sep=""),as.character(best_EC50))))
      pl_ci = as.data.frame(rbind(pl_ci,cbind(paste("slope",br,sep=""), as.character(best_slope))))
      createSheet(res.wb, name = paste(s,br,"_EC50",sep=""))
      createName(res.wb, name = paste(s,br,"_EC50",sep=""), formula = paste(s,br,"_EC50","!$B$2",sep = ""))
      writeWorksheet(res.wb,pl_ci,paste(s,br,"_EC50",sep=""))# save coefficients of each model (EC50)
      
      # save results into excel file (concentration values and viability for the fitted curves, from min to max, in 150 steps)
      cat(paste("Saving fitted curve to Workbook >>", b,"  Worksheet >>",s, "\n", sep=''))
      concentrations = seq(min(data_for_fit$conc),max(data_for_fit$conc),length.out=150)
      curve_values = predict(curvefit, newdata = data.frame(conc = c(concentrations)))
      curve_data = cbind(concentrations, curve_values)
      createSheet(res.wb, name = paste(s,br,"_Data",sep=""))
      createName(res.wb, name = paste(s,br,"_Data",sep=""), formula = paste(s,br,"_Data","!$B$2",sep = ""))
      writeWorksheet(res.wb,curve_data,paste(s,br,"_Data",sep=""))# save coefficients of each model (EC50)
      
      # make figure
      cat(paste("Saving figure to Workbook >>", b, "\n", sep=''))
      createSheet(res.wb, name = paste(s,br,"_drc",sep=""))
      createName(res.wb, name = paste(s,br,"_drc",sep=""), formula = paste(s,br,"_drc","!$B$2",sep = ""))
      png(file = paste("drc",".png",sep=""),width=1000,height=700)
      devAskNewPage(ask = FALSE)
      plot(data_for_fit$conc,data_for_fit$pViab, xlim= c(min(data_for_fit$conc),max(data_for_fit$conc)), ylim= c(-10,125))
      x_values = seq(min(data_for_fit$conc),max(data_for_fit$conc),length.out=1000)
      lines(c(x_values),predict(curvefit, newdata = data.frame(conc = c(x_values))))
      dev.off()
      addImage(res.wb,paste("drc",".png",sep=""),name=paste(s,br,"_drc",sep=""),originalSize=TRUE)
      
      # make profile likelihood figure
      cat(paste("Saving PL figure to Workbook >>", b, "\n", sep=''))
      createSheet(res.wb, name = paste(s,br,"_PL",sep=""))
      createName(res.wb, name = paste(s,br,"_PL",sep=""), formula = paste(s,br,"_PL","!$B$2",sep = ""))
      png(file = paste("figPL",".png",sep=""),width=1000,height=700)
      devAskNewPage(ask = FALSE)
      plot(conc_opt, goodnessfit)
      dev.off()
      addImage(res.wb,paste("figPL",".png",sep=""),name=paste(s,br,"_PL",sep=""),originalSize=TRUE)
      saveWorkbook(res.wb)
      
    }
  }
}
