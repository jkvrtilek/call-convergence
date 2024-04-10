# create food sharing rates from Gerry's PhD data

library(tidyverse)

# clear workspace
rm(list=ls())

# first download "vamp_donations19.csv" and "vamp_possible donors2.csv" from figshare

# set working directory to location of saved files
setwd("/Users/jkvrtilek/Desktop/OSU/PhD/GitHub/call-convergence/social_data/2015_cleaning")


# import donations
d <- read.csv("vamp_donations19.csv")
str(d)

# import possible donations
pd <- 
  read.csv("vamp_possible donors2.csv") %>% 
  # convert to matrix
  as.matrix() %>% 
  # convert to vector (stack the columns)
  as.vector() %>% 
  # back to df
  as_tibble() %>% 
  # mark the missing donors to remove later
  mutate(donor = ifelse(value=="", "DELETE", value)) %>% 
  mutate(rate =0) %>% 
  dplyr::select(-value)  

# check dimensions, nrow should be 36 donors x 183 donors lists = 6588
str(pd)

# get vector of one recipient for each donor list
t1 <- 
  d %>% 
  group_by(possible.donor.list, receiver) %>% 
  summarize(n=n()) %>% 
  pull(receiver)
# check its 183 donor lists
str(t1)  

# now repeat each entry 36 times
t2 <- rep(t1, each= 36)
str(t2) 

# combine recipients to get all possible donations
pd$recipient <- t2

# reorder columns
pd <- pd %>% dplyr::select(donor, receiver= recipient, rate)
str(pd)

# get actual donations: donor, recipient, and food-sharing rate 
ad <- 
  d %>% 
  dplyr::select(donor, receiver, rate= duration.rate) %>% 
  mutate(donor = ifelse(donor=="", "DELETE", donor)) %>% 
  as_tibble()

# combine possible and actual donations
rates2015 <- 
  rbind(pd, ad) %>% 
  group_by(donor, receiver) %>% 
  summarize(rate= mean(rate, na.rm=T), .groups= 'drop') %>% 
  filter(donor != "DELETE") %>% 
  filter(receiver != "DELETE") %>% 
  filter(donor != receiver)

# choose working directory for saving file
setwd("~/Dropbox/Dropbox/_working/_ACTIVE/__students/Julia Vrtilek/2022/kinship-call-similarity")

# save file
write.csv(rates2015, file= "2015_foodsharing_rates.csv")
  

