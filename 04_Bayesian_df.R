# combine data from 2015, 2017, 2019, and call distance to make a dataframe
# columns: Directed dyad, Vocal distance, Kinship, Allogrooming given, Allogrooming received, Food given, Food received, Capture population bat 1, Capture population bat 2, Sex 1, Sex 2
# Julia Vrtilek
# 29 March 2024, updated 10 April 2024, 10 May 2024

# set working directory
setwd("/Users/jkvrtilek/Desktop/OSU/PhD/GitHub/call-convergence/")

# load packages
library(tidyverse)
library(gdata)

# load kinship data
kin_raw <- read_csv("social_data/kinship-all-bats.csv") %>% 
  mutate(dir.dyad = paste(bat1,bat2,sep="_")) %>% 
  mutate(kin.r = kinship) %>% 
  select(dir.dyad,bat1,bat2,kin.r,5:8)

# change kinship for bats from same site (Tole) in diff seasons from 0 to NA
temp <- kin_raw %>% 
  filter(site1==site2) %>% 
  filter(season1!=season2)

kin <- kin_raw %>% 
  mutate(kinship = case_when(dir.dyad %in% temp$dir.dyad ~ NA,
                              !(dir.dyad %in% temp$dir.dyad) ~ kin.r)) %>% 
  select(dir.dyad:bat2, kinship, site1:season2)

# load phd data
rates15 <- read.csv("social_data/2010-2014_foodsharing/food_sharing_rates_2010-2014.csv") %>% 
  mutate(dir.dyad = paste(donor,receiver,sep="_")) %>% 
  mutate(grooming = NA) %>% 
  mutate(foodsharing = donation.rate) %>% 
  mutate(actor.sex = donor.sex) %>% 
  select(dir.dyad,grooming,foodsharing,actor.sex,receiver.sex)

# load and combine postdoc data
food17 <- read.csv("social_data/2016-2017_relationship_formation/vamp_dyadic_foodsharing_rates_2016-2017.csv") %>% 
  mutate(dir.dyad = paste(actor,receiver,sep="_")) %>% 
  mutate(foodsharing = donation.rate) %>% 
  select(dir.dyad,foodsharing)

groom17 <- read.csv("social_data/2016-2017_relationship_formation/vamp_dyadic_grooming_rates_2016-2017.csv") %>% 
  mutate(dir.dyad = paste(actor,receiver,sep="_")) %>% 
  mutate(grooming = grooming.rate) %>% 
  select(dir.dyad,grooming)

rates17 <- left_join(groom17,food17,by="dir.dyad")

# load lab data
rates19 <- read.csv("social_data/2019_grooming/2019_grooming_rates.csv") %>% 
  mutate(grooming = grooming.rate) %>% 
  mutate(foodsharing = NA) %>% 
  select(dir.dyad,grooming,foodsharing,actor.sex,receiver.sex)

# combine rates for left_join
all_rates <- rbind(rates15,rates17,rates19)

# load DFA distances
dist <- as.matrix(read_csv("vocal-distance-lda.csv"))
rownames(dist) <- colnames(dist)
distdf <- enframe(unmatrix(dist))
colnames(distdf) <- c("dir.dyad","dist")
distdf$dir.dyad <- gsub(":","_",distdf$dir.dyad)

# make df with all extant data
d <- left_join(kin,all_rates,by = "dir.dyad")
d2 <- left_join(d,distdf,by = "dir.dyad")

d3 <- d2 %>% 
  separate(dir.dyad, into = c(actor,receiver), sep = "_", remove = FALSE)

write.csv(d2,"vocal_social_data.csv")
