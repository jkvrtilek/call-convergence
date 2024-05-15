# combine data from 2015, 2017, 2019, and call distance to make a dataframe
# columns: Directed dyad, Vocal distance, Kinship, Allogrooming, Foodsharing, Capture population bat 1, Capture population bat 2, Sex 1, Sex 2
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
  select(dir.dyad, kinship, site1:season2)

# load phd data
rates15 <- read.csv("social_data/2010-2014_foodsharing/food_sharing_rates_2010-2014.csv") %>% 
  mutate(dir.dyad = paste(donor,receiver,sep="_")) %>% 
  mutate(grooming = NA) %>% 
  mutate(foodsharing = donation.rate) %>%
  mutate(actor.sex = toupper(actor.sex)) %>% 
  mutate(receiver.sex = toupper(receiver.sex)) %>% 
  select(dir.dyad,foodsharing,grooming,actor.sex,receiver.sex)

# load and combine postdoc data
food17 <- read.csv("social_data/2016-2017_relationship_formation/vamp_dyadic_foodsharing_rates_2016-2017.csv") %>% 
  mutate(actor = case_when(actor == "son.of.ola" ~ "son-of-ola",
                           actor == "no.band" ~ "no-band",
                           .default = actor)) %>% 
  mutate(receiver = case_when(receiver == "son.of.ola" ~ "son-of-ola",
                              receiver == "no.band" ~ "no-band",
                              .default = receiver)) %>% 
  mutate(dir.dyad = paste(actor,receiver,sep="_")) %>% 
  mutate(foodsharing = donation.rate) %>% 
  select(dir.dyad,foodsharing)

groom17 <- read.csv("social_data/2016-2017_relationship_formation/vamp_dyadic_grooming_rates_2016-2017.csv") %>% 
  mutate(actor = case_when(actor == "son.of.ola" ~ "son-of-ola",
                           actor == "no.band" ~ "no-band",
                           .default = actor)) %>% 
  mutate(receiver = case_when(receiver == "son.of.ola" ~ "son-of-ola",
                           receiver == "no.band" ~ "no-band",
                           .default = receiver)) %>% 
  mutate(dir.dyad = paste(actor,receiver,sep="_")) %>% 
  mutate(grooming = grooming.rate) %>% 
  select(dir.dyad,grooming,actor.sex,receiver.sex)

rates17 <- left_join(food17,groom17,by="dir.dyad")

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

# remove bat w/self
distdf2 <- distdf %>% 
  separate(dir.dyad, into = c("bat1","bat2"), sep = "_", remove = FALSE) %>% 
  filter(bat1 != bat2) %>% 
  select(-bat1,-bat2)

# make df with all extant data
d <- left_join(kin,all_rates,by = "dir.dyad")
d2 <- left_join(d,distdf2,by = "dir.dyad") %>% 
  filter(!is.na(dist))
d3 <- d2 %>% 
  separate(dir.dyad, into = c("actor","receiver"), sep = "_", remove = FALSE)

# fill in sexes
m <- c("vampison","black-panther","falcon","hulk","bat5","jenna-bat")
f <- c("bat6","bat7","lds","scc","six")

sex <- d3 %>% 
  select(actor, actor.sex) %>% 
  distinct() %>% 
  mutate(sex = case_when(actor %in% m ~ "M",
                         actor %in% f ~ "M",
                         .default = actor.sex)) %>% 
  select(actor, sex) %>% 
  filter(!is.na(sex))

female <- sex %>% 
  filter(sex == "F")

male <- sex %>% 
  filter(sex == "M")

d4 <- d3 %>% 
  mutate(a.sex = case_when(actor %in% female$actor ~ "F",
                           actor %in% male$actor ~ "M",
                           .default = NA)) %>% 
  mutate(r.sex = case_when(receiver %in% female$actor ~ "F",
                           receiver %in% male$actor ~ "M",
                           .default = NA)) %>% 
  select(dir.dyad:receiver, a.sex, r.sex, site1:season2, kinship, grooming, foodsharing, dist)

write.csv(d4,"vocal_social_data.csv")
