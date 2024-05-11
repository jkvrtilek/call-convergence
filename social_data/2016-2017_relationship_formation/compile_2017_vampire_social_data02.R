# compile all vampire bat food sharing data from Gerry Carter's postdoc

# download PhD data from here: 
# Carter, Gerald (2019). Data and R Code for Carter et al. 2020. Development of new food-sharing relationships in vampire bats. Current Biology.. figshare. Dataset. https://doi.org/10.6084/m9.figshare.11369268.v1

# edited by Julia, 8 May 2024

# clear workspace
rm(list=ls())

# load packages
library(tidyverse)

# set working directory
setwd("/Users/jkvrtilek/Desktop/OSU/PhD/GitHub/call-convergence/social_data/2016-2017_relationship_formation")

# get data
load("new_bonds_data.RData")

unique(rates$date)
unique(rates$days.together)

# get relevant data
d <- 
  rates %>% 
  as_tibble() %>% 
  # exclude controlled introductions?
  # filter(cage.type=="colony") %>% 
  filter(behav== "f"| behav== "g") %>%
  filter(!nonfocal) %>% 
  # clarify column meanings
  select(datetime, period, trial, trial.duration, subject, actor, receiver, behav, duration.sec= rate, actor.origin= a.pop, receiver.origin= r.pop, actor.is.mom, receiver.is.mom, new.dyad, both.adult= adult.dyad, intro.date, days.together, actor.birthdate= a.birth, receiver.birthdate= r.birth, kinship= kinship2)

# save daily rates
write.csv(d, "vamp_daily_rates_2016-2017.csv")

# get sexes
sex <- read.csv("pup_attributes.csv") 

# d11 is ddd, replace
sex$name <- gsub("d11","ddd",sex$name)

# make actor and receiver dfs to join
actor.sex <- sex %>% 
  mutate(actor = name) %>% 
  mutate(actor.sex = toupper(sex)) %>% 
  select(actor, actor.sex)

receiver.sex <- sex %>% 
  mutate(receiver = name) %>% 
  mutate(receiver.sex = toupper(sex)) %>% 
  select(receiver, receiver.sex)

# get dyadic grooming
g <- 
  d %>% 
  filter(behav=="g") %>% 
  group_by(actor, receiver) %>% 
  summarize(grooming.sec= sum(duration.sec, na.rm=T), observation.hours= sum(trial.duration)) %>% 
  ungroup() %>% 
  mutate(observation.seconds = observation.hours*60*60) %>% 
  mutate(grooming.rate = grooming.sec/observation.seconds)

g2 <- g %>% 
  left_join(actor.sex, by = "actor") %>% 
  left_join(receiver.sex, by = "receiver")
g2$actor.sex <- g2$actor.sex %>% replace_na("F")
g2$receiver.sex <- g2$receiver.sex %>% replace_na("F")

write.csv(g2, "vamp_dyadic_grooming_rates_2016-2017.csv")
  
# get dyadic food-sharing
f <- 
  d %>% 
  filter(behav=="f") %>% 
  group_by(actor, receiver) %>% 
  summarize(sharing.sec= sum(duration.sec, na.rm=T), observation.hours= sum(trial.duration)) %>% 
  ungroup() %>% 
  mutate(observation.seconds = observation.hours*60*60) %>% 
  mutate(donation.rate = sharing.sec/observation.seconds)

f2 <- f %>% 
  left_join(actor.sex, by = "actor") %>% 
  left_join(receiver.sex, by = "receiver")
f2$actor.sex <- f2$actor.sex %>% replace_na("F")
f2$receiver.sex <- f2$receiver.sex %>% replace_na("F")

write.csv(f2, "vamp_dyadic_foodsharing_rates_2016-2017.csv")

