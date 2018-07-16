# DATA DESCRIPTION
# In the ring study, 6 different laboratories have measured the dose-response curve in RTgill cell line exposed to 6 different chemicals, 
# using different reporter dyes. For each chemical they have measured the EC50 3-times (or 2 or 4). 
# It is important to evaluate how close their measurements are to see how much we can trust the essays.
#
# WHAT IS AN OUTLIER?
# an outlier is a data point (or a set of data points) that is very different from other datapoints in a dataset, so different in fact, 
# that it is believed to be the results of an error and needs to be removed from the dataset before the analysis is run. 

# LABORATORIES AS OUTLIERS
# We define laboratories as outliers when the measured mean EC50 of the individual laboratory is beyond X (THRESHOLD) times the standard deviation 
# among all laboratories performing the same experiment.

# BIOLOGICAL REPLICATES AS OUTLIER
# We define a biological replicate within a laboratory as an outlier if it's more than X (THRESHOLD) times the standard deviation of all 
# laboratories away from the laboratory mean 



# clear space
rm(list = ls())

# load libraries
library("reshape2")
library("readxl")


## Read the input file for outlier analysis
# path
filePathXXX_xlsx = "//eawag/userdata/zupanian/My Documents/anze/scripts/R scripts/outliers_ring/outlier_inputFile.xlsx"
# read file
data<- read_excel(filePathXXX_xlsx)

## ENTER THRESHOLD HERE
threshold = 2.0

## LABORATORY OUTLIERS
# format data (ignore column one where laboratory names are written)
EC50s = data[1:nrow(data),2:ncol(data)]  
# average EC50s by lab (row)
EC50s_mean = rowMeans(EC50s, na.rm = TRUE)
# calculate SD across all data
stdevi = sd(EC50s_mean)
# is there any laboratory EC50 beyond the threshold*sd?
laboratory_outlier_index = which(abs(EC50s_mean-mean(EC50s_mean))>threshold*stdevi)
# laboratory_outlier_index gives you the row index of the laboratory that is an outlier 


# BIOLOGICAL REPLICATE OUTLIERS
# format data (ignore column one where laboratory names are written)
EC50s = data[1:nrow(data),2:ncol(data)]  
# normalize data by mean EC50 for each laboratory
EC50s_norm = EC50s /rowMeans(EC50s, na.rm = TRUE, dims = 1)
# calculate SD across all data
vector_norm = melt(EC50s_norm, na.rm = TRUE)
stdevi = sd(vector_norm[,2])
# is there any biological replicates beyond the threshold*sd?
biological_outlier_index = which(abs(EC50s_norm-1)>threshold*stdevi, arr.ind = T)
# biological_outlier_index gives you the row and column index of the outliers 







