# extract 2019 grooming data from Imran's stuff
# Julia Vrtilek
# 29 March 2024, updated 7 April 2024 for published data

setwd("/Users/jkvrtilek/Desktop/OSU/PhD/GitHub/call-convergence/social_data/")

library(tidyverse)
library(gdata)

# get grooming data
load("2019_cleaning/rates2019_2021-05-18.Rdata")

rates <- rates2019 %>% 
  filter(cage == "big_cage") %>% 
  filter(behav == "g") %>% 
  select(6:11)

# get key
key <- read_csv("2019_cleaning/2019batkey.csv")

# change actor name
key_act <- key %>% 
  mutate(actor = final.bat.ID) %>%
  mutate(actor2 = bat.name) %>% 
  select(actor,actor2)

rates2 <- left_join(rates, key_act, by = "actor")

# change bat2 name
key_rec <- key %>% 
  mutate(receiver = final.bat.ID) %>%
  mutate(receiver2 = bat.name) %>% 
  select(receiver,receiver2)

rates3 <- left_join(rates2, key_rec, by = "receiver")

rates4 <- rates3 %>% 
  mutate(dir.dyad = paste(actor2,receiver2,sep="_")) %>% 
  select(dir.dyad,actor2,receiver2,behav,rate) %>% 
  group_by(dir.dyad, behav) %>% 
  summarize(rate= mean(rate, na.rm=T))

saveRDS(rates4,"2019_grooming_rates.RDS")
