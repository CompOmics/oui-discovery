#!/usr/bin/env

# install.packages("devtools")
# devtools::install_github("freejstone/groupwalk")

# parse arguments
args <- commandArgs(trailingOnly = TRUE)
dataset_dir <- args[1]
working_dir <- args[2]
file_name <- args[3]

##@
##working_dir <- "D:/vasv1101/Documents_discC/tpp_pride_reanalysis"
##dataset_dir <- paste0(working_dir,"/tpp_runs/test_groupwalk")

# set working directory
setwd(working_dir) 
# in RStudio, can be set with Session -> set working directory -> set to source file location

# read data file
##data <- read.csv(paste(dataset_dir,"/group-walk-input.csv", sep=''))
#@
data <- read.csv(paste(dataset_dir,"/", file_name , sep=''))
data$FDRGroup <- as.factor(data$FDRGroup)
data$isTarget <- data$isTarget == "True"

# calculate FDR
source("Group_walk.R")
data$group_qval_fval <- group_walk(data$fval, data$isTarget, data$FDRGroup)
##data$group_qval_pep <- group_walk(data$PEP, data$isTarget, data$FDRGroup)
data$group_qval_prob <- group_walk(data$peptideprophet_probability, data$isTarget, data$FDRGroup)
##data$group_qval_xcorr <- group_walk(data$xcorr, data$isTarget, data$FDRGroup)
write.table(data, paste(dataset_dir,"/groupwalk_output_",file_name, sep=''), row.names = FALSE, sep=",")
