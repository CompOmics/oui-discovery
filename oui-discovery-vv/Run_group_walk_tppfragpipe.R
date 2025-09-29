#!/usr/bin/env

# install.packages("devtools")
# devtools::install_github("freejstone/groupwalk")

# parse arguments
args <- commandArgs(trailingOnly = TRUE)
dataset_dir <- args[1]
working_dir <- args[2]
file_name <- args[3]

# set working directory
setwd(working_dir)

# read data file
data <- read.csv(paste(dataset_dir,"/", file_name , sep=''))
data$FDRGroup <- as.factor(data$FDRGroup)
data$isTarget <- data$isTarget == "True"

# calculate FDR
source("Group_walk.R")
data$group_q_prob <- group_walk(data$peptideprophet_probability, data$isTarget, data$FDRGroup)
write.table(data, paste(dataset_dir,"/groupwalk_output_",file_name, sep=''), row.names = FALSE, sep=",")
