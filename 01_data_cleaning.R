# data cleaning: relabel files, s to ms, duration and peak freq filters, add binary ff variable
# Julia Vrtilek
# 23 March 2023

# set working directory
setwd("/Users/jkvrtilek/Desktop/OSU/PhD/GitHub/call-convergence")

# load packages
library(tidyverse)

# load data
raw <- read.csv("vampire_call_measures.csv")

# rename files mistakenly labeled from 2018
d <- raw
length(grep("2018-", raw$sound.files))
d$sound.files <- gsub("2018-08-29","2017-09-03", raw$sound.files)
length(grep("2018-", d$sound.files))

# convert seconds to milliseconds
d2 <- d
d2$duration <- d$duration * 1000

# duration and peak frequency filters
d3 <- d2 %>% 
  filter(duration > 3) %>% 
  filter(duration < 50) %>% 
  filter(peakf > 10) %>% 
  filter(peakf < 30)

# add binary variable for whether fundamental frequency measures succeeded
d4 <- d3 %>% 
  mutate(indicator = case_when(!is.na(meanslope) ~ T,
                               is.na(meanslope) ~ F,
                               is.nan(meanslope) ~ F))

# save file
write.csv(d4, "vampire_call_measures_filtered.csv")
