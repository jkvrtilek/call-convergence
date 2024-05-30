# compile all vampire bat food sharing data from Gerry Carter's PhD
# Gerry Carter and Julia Vrtilek

# download Michigan/OBC data from here: 
# Carter, Gerald; Schino, Gabriele; Farine, Damien (2018). Data and R code for "Challenges in assessing the roles of nepotism and reciprocity in cooperation networks". figshare. Dataset. https://doi.org/10.6084/m9.figshare.6072272.v6 

# clear workspace
rm(list=ls())

# load packages
library(tidyverse)

# set working directory
setwd("/Users/jkvrtilek/Desktop/OSU/PhD/GitHub/call-convergence/social_data/2010-2014_foodsharing")


# wrangle and sanity check UMD colony data -----

# get data from University of Maryland (UMD) colony (from unpublished oxytocin pilot study)
umd <- 
  read.csv("2013_oxytocin_social_data.csv") %>% 
  mutate(subject= paste0("bat",subject)) %>% 
  mutate(donor= paste0("bat",donor)) %>% 
  # numbers are count of minutes that have at least 5 sec of behavior (multiplying by 60 will overestimate actual seconds)
  # convert minutes to seconds assuming 40 s of sharing per minute with at least 5 s of food sharing
  # this adjustment creates matching food transfer rates based on all other data (see linear model below)
  mutate(grooming = (grooming*40), sharing = (sharing*40)) %>% 
  as_tibble()

# check that sharing per trial causes mass gain to a degree that matches past work
t <- 
  umd %>% 
  group_by(trial) %>% 
  summarize(mass.gain= mean(trial.mass.gain.subject, na.rm=T), sharing= sum(sharing, na.rm=T)) %>% 
  mutate(sharing.min= sharing/60, mass.gain.mg= mass.gain*1000)

# plot looks similar to other work (see Figure S1 in Carter & Wilkinson 2015 Proc B)
t %>%   
  ggplot(aes(x=sharing, y=mass.gain))+
    geom_point(size=2)+
    geom_smooth(method= "lm")+
  xlab("mouthlicking (seconds)")+
  ylab("weight gain (grams)")

# expected value from my other PhD data (121 fasting trials) is that 37 milligrams of food [95% CI: 31–43 mg] should be transferred per minute of mouth licking
lm(mass.gain.mg~sharing.min, data=t)
# observed value is 37 mg per minute


# wrangle and sanity check Michigan/OBC data -----

# get data from Michigan colony
# food sharing events
d <- read.csv("vamp_donations19.csv") 
# possible donors for each trial
pos.d <- read.csv("vamp_possible donors3.csv") 
# kinship estimates using two methods (with and without pedigree info)
kinship <- read.csv("vamp_relatedness12.csv") %>% 
  mutate(dyad = paste(bat1, bat2))
# bat sexes and ages
bats <- read.csv("vamp_sex_age.csv")

# check all receivers and donors are in bats dataframe (df)
unique(d$receiver)[which(!unique(d$receiver) %in% bats$name)]
unique(d$donor)[which(!unique(d$donor) %in% bats$name)]
# fix bat names
d$donor[which(d$donor== "leord")] <- "leonard"
d$donor[which(d$donor== "")] <- NA

# check which bats are not in Michigan donation data
tbats <- unique(c(d$receiver, d$donor))
bats$name[which(!bats$name %in% tbats)]
# these are bats that never donated in Michigan

# get names in possible donors that don't match bat names
pd <- unique(as.vector(as.matrix(pos.d)))
pd[which(!pd %in% bats$name)]

# change these names in possible donors matrix
pos.d2 <- pos.d %>% 
  mutate(across(everything(), ~ str_replace(string = .x, pattern = "two", replacement = "bat2"))) %>% 
  mutate(across(everything(), ~ str_replace(string = .x, pattern = "three", replacement = "bat3"))) %>% 
  mutate(across(everything(), ~ str_replace(string = .x, pattern = "four", replacement = "bat4")))

# update possible donors
pd <- unique(as.vector(as.matrix(pos.d2)))
bats$name[which(!bats$name %in% pd)]
# these are bats that were never tested as donors in Michigan

# to get actual and possible donations...
# first get all combinations of donors and trials
t <- 
  expand_grid(trial_receiver= unique(paste(d$trial, d$receiver)), donor= unique(bats$name), duration= 0) %>% 
  separate(trial_receiver, into= c("trial", "receiver")) %>% 
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

# save food-sharing trial data (Michigan/OBC colony)  
# write.csv(d2, "food_donations_OBC_2010-2014.csv")

# combine Michigan/OBC and UMD data into dyadic df -----

# get dyadic food sharing rates from Michigan colony
d3 <- 
  d2 %>% 
  group_by(donor, receiver) %>% 
  summarize(donation.sec= sum(donation.sec), possible.sec= sum(observation.sec), .groups= 'drop') 

# get dyadic food sharing rates from UMD colony
umd2 <- 
  umd %>% 
  mutate(receiver= subject) %>% 
  mutate(possible.sec= 60*60) %>% 
  group_by(donor, receiver) %>% 
  summarize(
    possible.sec= sum(possible.sec, na.rm=T), 
    donation.sec = sum(sharing, na.rm=T), 
    .groups = 'drop')
  
# combine and add other dyadic variables
d4 <- 
  d3 %>% 
  full_join(umd2) %>% 
  group_by(donor, receiver) %>% 
  summarize(
    possible.sec= sum(possible.sec, na.rm=T), 
    donation.sec = sum(donation.sec, na.rm=T), 
    .groups = 'drop') %>% 
  mutate(dyad= paste(donor,receiver)) %>%  
  # add kinship
  mutate(kinship= kinship$pedigree.m19[match(.$dyad, kinship$dyad)]) %>% 
  # add donor sex
  mutate(donor.sex= bats$sex[match(.$donor, bats$name)]) %>% 
  # add receiver sex
  mutate(receiver.sex= bats$sex[match(.$receiver, bats$name)]) %>% 
  # add donor age
  mutate(donor.age.2014.10.07= bats$age20141007[match(.$donor, bats$name)]) %>% 
  # add receiver age
  mutate(receiver.age.2014.10.07= bats$age20141007[match(.$receiver, bats$name)]) %>% 
  # donation rate
  mutate(donation.rate= donation.sec/possible.sec) %>% 
  select(-dyad) 

# save food-sharing rates data  
write.csv(d4, "food_sharing_rates_2010-2014.csv")
