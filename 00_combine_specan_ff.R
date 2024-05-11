# script to combine the spectro_analysis results and the fundamental frequency summary measures
# Julia Vrtilek
# 21 March 2024

# set working directory
setwd("/Users/jkvrtilek/Desktop/OSU/PhD/GitHub/call-convergence")

# load packages
library(tidyverse)

# read in spectro_analysis data
specan <- readRDS("ds_completespecan.RDS")

# read in fund freq data
ff <- readRDS("ds_fundfreq_summary.RDS")

# combine datasets by sound.files
fulld <- left_join(specan, ff, by = "sound.files")

# save combined dataset
write.csv(fulld, "vampire_call_measures.csv")
