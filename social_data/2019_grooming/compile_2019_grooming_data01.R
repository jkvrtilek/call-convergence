# extract 2019 grooming data from Imran's stuff
# Julia Vrtilek
# revised 9 May 2024

# clear workspace
rm(list=ls())

# load packages
library(tidyverse)
library(gdata)

setwd("/Users/jkvrtilek/Desktop/OSU/PhD/GitHub/call-convergence/social_data/2019_grooming")

# get grooming data
load("rates2019_2021-05-18.Rdata")

rates <- rates2019 %>% 
  filter(cage == "big_cage") %>% 
  filter(behav == "g") %>% 
  select(6:11) %>% 
  group_by(dyad,actor,receiver) %>% 
  summarize(grooming.sec = sum(rate), hours = n()) %>% 
  ungroup() %>% 
  mutate(seconds = hours*60*60) %>% 
  mutate(grooming.rate = grooming.sec/seconds)

# get key
key <- read_csv("2019batkey.csv")

# change actor name
key_act <- key %>% 
  mutate(actor = final.bat.ID) %>%
  mutate(actor2 = bat.name) %>% 
  mutate(actor.sex = sex) %>% 
  select(actor,actor2,actor.sex)

rates2 <- left_join(rates, key_act, by = "actor")

# change bat2 name
key_rec <- key %>% 
  mutate(receiver = final.bat.ID) %>%
  mutate(receiver2 = bat.name) %>% 
  mutate(receiver.sex = sex) %>% 
  select(receiver,receiver2,receiver.sex)

rates3 <- left_join(rates2, key_rec, by = "receiver")

rates4 <- rates3 %>% 
  mutate(dir.dyad = paste(actor2,receiver2,sep="_")) %>% 
  select(dir.dyad,actor2,receiver2,grooming.rate,actor.sex,receiver.sex) %>% 
  mutate(actor = actor2) %>% 
  mutate(receiver = receiver2) %>% 
  select(dir.dyad,actor,receiver,grooming.rate,actor.sex,receiver.sex)

write.csv(rates4,"2019_grooming_rates.csv")
