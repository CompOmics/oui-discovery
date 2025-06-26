#!/usr/bin/env

# install.packages("devtools")
# devtools::install_github("freejstone/groupwalk")

# parse arguments
args <- commandArgs(trailingOnly = TRUE)
working_dir <- args[1]
input       <- args[2]
output      <- args[3]

# set working directory
setwd(working_dir) 
# in RStudio, can be set with Session -> set working directory -> set to source file location

# read data file
data <- read.csv(input)
data$FDRGroup <- as.factor(data$FDRGroup)

# calculate FDR
source("Group_walk.R")
data$group_qval <- group_walk(data$psm_score, data$isTarget, data$FDRGroup)
write.table(data, output, row.names = FALSE, sep=",")
