# compile all vampire bat food sharing data from Gerry Carter's PhD

# download PhD data from here: 
# Carter, Gerald; Schino, Gabriele; Farine, Damien (2018). Data and R code for "Challenges in assessing the roles of nepotism and reciprocity in cooperation networks". figshare. Dataset. https://doi.org/10.6084/m9.figshare.6072272.v6 

# edited by Julia, 8 May 2024

rm(list=ls())

# load packages
library(tidyverse)

# set working directory
setwd("/Users/jkvrtilek/Desktop/OSU/PhD/GitHub/call-convergence/social_data/2010-2014_foodsharing")

# get data
# food sharing events
d <- read.csv("vamp_donations19.csv") 
# possible donors for each trial
pos.d <- read.csv("vamp_possible donors3.csv") 
# kinship estimates using two methods (with and without pedigree info)
kinship <- read.csv("vamp_relatedness12.csv") 
# bat sexes and ages
bats <- read.csv("vamp_sex_age.csv")

# check all receivers and donors are in bats dataframe (df)
unique(d$receiver)[which(!unique(d$receiver) %in% bats$name)]
unique(d$donor)[which(!unique(d$donor) %in% bats$name)]
# fix bat names
d$donor[which(d$donor== "leord")] <- "leonard"
d$donor[which(d$donor== "")] <- NA

# check which bats are not in donation data
tbats <- unique(c(d$receiver, d$donor))
bats$name[which(!bats$name %in% tbats)]
# correct: these are all bats that just never donated

# get names in possible donors that don't match bat names
pd <- unique(as.vector(as.matrix(pos.d)))
pd[which(!pd %in% bats$name)]

# change these names in possible donors matrix
pos.d2 <- pos.d %>% 
  mutate(across(everything(), ~ str_replace(string = .x, pattern = "two", replacement = "bat2"))) %>% 
  mutate(across(everything(), ~ str_replace(string = .x, pattern = "three", replacement = "bat3"))) %>% 
  mutate(across(everything(), ~ str_replace(string = .x, pattern = "four", replacement = "bat4")))

bats$name[which(!bats$name %in% pd)]
# correct: these are all bats that were never tested as donors

# to get actual and possible donations...
# first get all combinations of bats for every trial
t <- 
  expand_grid(trial= d$trial, donor= bats$name, receiver= bats$name, duration= 0) %>% 
  filter(donor!=receiver)

# then delete impossible cases of food sharing
for (i in 1:length(pos.d2)) {
  # get current trial
  focal.trial <- i
  # get possible donors for current trial (bats present during trial)
  pos <- pos.d2[,focal.trial]
  # get the possible receiver (fasted bat) for current trial
  pos2 <- 
    d %>% 
    filter(trial== focal.trial) %>% 
    pull(receiver) %>% 
    unique()
  # put NA duration for trials where donor is not a possible donor or receiver is not possible receiver
  t$duration[which(t$trial == focal.trial & !(t$donor %in% pos))] <- NA
  t$duration[which(t$trial == focal.trial & !(t$receiver %in% pos2))] <- NA
  
  # show progress
  print(paste(i,"of",length(pos.d2)))
}

# label observed cases of food sharing
d$trial.donor.receiver <- paste0(d$trial, d$donor, d$receiver)

# get actual and possible donations
d2 <- 
  t %>% 
  # discard impossible donations
  filter(duration==0) %>% 
  # add observed cases of food sharing
  mutate(trial.donor.receiver= paste0(trial,donor,receiver)) %>% 
  mutate(donation.sec = d$duration[match(.$trial.donor.receiver, d$trial.donor.receiver)]) %>% 
  mutate(donation.sec= ifelse(is.na(donation.sec), 0, donation.sec)) %>% 
  # add trial duration (observation time) in seconds
  mutate(observation.sec = 60*60*d$trial.duration..h.[match(.$trial, d$trial)]) %>% 
  # add trial dates
  mutate(date = d$date[match(.$trial, d$trial)]) %>% 
  select(date, trial, donor, receiver, donation.sec, observation.sec)

# save food-sharing trial data  
write.csv(d2, "food_sharing_donations_2010-2014.csv")

# label kinship dyads
kinship$dyad <- paste(kinship$bat1, kinship$bat2)

# now get dyadic rates with sex, kinship
d3 <- 
  d2 %>% 
  group_by(donor, receiver) %>% 
  summarize(donation.sec= sum(donation.sec), possible.sec= sum(observation.sec), .groups= 'drop') %>% 
  # label dyad
  mutate(dyad= paste(donor,receiver)) %>% 
  # add kinship
  mutate(kinship= kinship$pedigree.m19[match(.$dyad, kinship$dyad)]) %>% 
  # add donor sex
  mutate(actor.sex= bats$sex[match(.$donor, bats$name)]) %>% 
  # add receiver sex
  mutate(receiver.sex= bats$sex[match(.$receiver, bats$name)]) %>% 
  # donation rate
  mutate(donation.rate= donation.sec/possible.sec) %>% 
  select(-dyad) 

# save food-sharing rates data  
write.csv(d3, "food_sharing_rates_2010-2014.csv")
