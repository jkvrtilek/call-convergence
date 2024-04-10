# re-arrange data published in Current Biology 2020
# Gerry Carter
# 6 April 2024

# set working directory
setwd("/Users/jkvrtilek/Desktop/OSU/PhD/GitHub/call-convergence")

library(tidyverse)
library(igraph)
library(gdata)

# function to convert a list of dyadic interactions to a network ------
a_b_edgelist_to_matrix <- function(el=el, symbol="_", directed= T, make.NA.zero=T){
  a <- str_split(as.data.frame(el)[,1],symbol, simplify = TRUE)[,1]
  r <- str_split(as.data.frame(el)[,1],symbol, simplify = TRUE)[,2]
  y <- as.data.frame(el)[,2]
  e <- data.frame(a,r,y, stringsAsFactors = F)
  require(igraph)
  if (make.NA.zero){
    g <- graph_from_data_frame(e, directed=directed)
    m <- get.adjacency(g, attr='y', sparse=FALSE)
    m
  }else{
    e$y <- e$y+1 # temporarily add one to distinguish between 0 and NA
    g <- graph_from_data_frame(e, directed=directed)
    m <- get.adjacency(g, attr='y', sparse=FALSE)
    m[m==0] <- NA # relabel missing values as NA
    m <- m-1 # subtract one to adjust values back
    m
  }
}

# get data
load("social_data/2017_cleaning/new_bonds_data.RData")

# fix names
rates$dir.dyad <- gsub("no.band","no-band",rates$dir.dyad)
rates$dir.dyad <- gsub("son.of.ola","son-of-ola",rates$dir.dyad)


# get food sharing network
feed.net <- 
  rates %>% 
  filter(period == "big cage") %>% 
  filter(behav== 'f') %>% 
  group_by(dir.dyad, behav) %>% 
  summarize(rate= mean(rate, na.rm=T))

saveRDS(feed.net,"social_data/2017_foodsharing_rates.RDS")


# get allogrooming network
groom.net <- 
  rates %>% 
  filter(period == "big cage") %>% 
  filter(behav== 'g') %>% 
  group_by(dir.dyad, behav) %>% 
  summarize(rate= mean(rate, na.rm=T)) 

saveRDS(groom.net,"social_data/2017_grooming_rates.RDS")