# setup ------

# clear work space
rm(list=ls())

# load packages
library(tidyverse)
library(boot)
library(asnipe)
library(vegan)
library(lme4)
library(lmerTest)
library(glmmTMB)
library(broom.mixed)
library(brms)
library(performance)
library(igraph)
library(bayesplot)

# set working directory
setwd("/Users/jkvrtilek/Desktop/OSU/PhD/GitHub/call-convergence")

# get data
d <- read.csv("vocal_social_data.csv") %>% 
  select(-X)

# add undirected dyads
d2 <- d %>% 
  mutate(dyad = ifelse(actor<receiver, paste(actor,receiver, sep= "_"), paste(receiver,actor, sep= "_"))) %>% 
  arrange(dyad) %>% 
  group_by(dyad) %>% 
  summarize(kinship = mean(kinship, na.rm=T),
            grooming = mean(grooming, na.rm=T),
            sharing = mean(foodsharing, na.rm=T),
            distance = mean(dist, na.rm=T)) %>% 
  separate(dyad, into= c("actor", "receiver"), sep= "_", remove=F)

# restore site and season
seasonsite <- d %>% 
  mutate(dyad = dir.dyad) %>% 
  select(dyad,a.sex,r.sex,site1,site2,season1,season2)

d3 <- d2 %>% 
  left_join(seasonsite, by = "dyad") %>% 
  # add other social dyadic variables
  # add capture site
  mutate(site= case_when(
    site1 == "zoo" & site2 == "zoo" ~ "US",
    site1 == "chorrera" & site2 == "chorrera" ~ "CH",
    site1 == "tole" & site2 == "tole" ~ "TL",
    site1 == "las.pavas" & site2 == "las.pavas" ~ "LP",
    site1 == "lake.bayano" & site2 == "lake.bayano" ~ "LB",
    site1 != site2 ~ "different sites",
    TRUE ~ "other")) %>% 
  # add lab colony
  mutate(season= case_when(
    season1 == "phd" & season2 == "phd" ~ "2015",
    season1 == "postdoc" & season2 == "postdoc" ~ "2017",
    season1 == "lab" & season2 == "lab" ~ "2019",
    season1 != season2  ~ "different colonies",
    TRUE ~ "other")) %>% 
  # add zoo vs wild
  mutate(site.type= case_when(
    site1 == "zoo" & site2 == "zoo" ~ "zoo",
    site1 == "zoo" & site2 != "zoo" ~ "different",
    site1 != "zoo" & site2 == "zoo" ~ "different",
    TRUE ~ "wild")) %>% 
  mutate(season_site= paste(season, site, sep="_")) %>% 
  mutate(treatment= case_when(
    season_site == "2019_CH" ~ "same wild &\n same lab colony",
    season_site == "2019_LB" ~ "same wild &\n same lab colony",
    season_site == "2017_LP" ~ "same wild &\n same lab colony",
    season_site == "2019_TL" ~ "same wild &\n same lab colony",
    season_site == "2017_TL" ~ "same wild &\n same lab colony",
    season_site == "2015_US" ~ "same long-term lab colony",
    season_site == "2017_different sites" ~ "different wild colony &\n same lab colony",
    season_site == "2019_different sites" ~ "different wild colony &\n same lab colony",
    season_site == "different colonies_TL" ~ "same wild colony &\n different lab colony",
    season_site == "different colonies_different sites" ~ "different wild colony &\n different lab colony",
    TRUE ~ "other"))

# functions ------

# convert a list of dyadic interactions to a network
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

# foodsharing ------
foodsharing <- 
  d3 %>% 
  filter(!is.na(distance)) %>% 
  filter(!is.na(sharing)) %>% 
  filter(kinship < 0.05) %>% # non-kin only
  filter(a.sex == "F") %>% 
  filter(r.sex == "F")


# model dist as predicted by food sharing: PHD ------
fs_phd <- foodsharing %>% 
  filter(season1=="phd" & season2== "phd")

# MRQAP
fsm1 <- 
  fs_phd %>% 
  select(dyad, sharing) %>% 
  a_b_edgelist_to_matrix(directed=F )
fdm1 <- 
  fs_phd %>% 
  select(dyad, distance) %>% 
  a_b_edgelist_to_matrix(directed=F)
mrqap.dsp(fdm1~fsm1)  
  
# choose number of chains and chain length
# adjust as necessary
nchains <- 4
chain_length <- 3000
warmup_length <- 1000

# fit brms model of effect of familiarity on vocal distance
fit_fs1 <-
  brm(scale(distance) ~ 
        scale(sharing) + 
        (1|mm(actor,receiver)),
      data = fs_phd, 
      family = "gaussian",
      cores = nchains,
      chains = nchains,
      iter = chain_length,
      warmup = warmup_length)

# check model
hist(residuals(fit_fs1))
pp_check(fit_fs1, ndraws=100)
summary(fit_fs1)

# get model results
coefs_fs1 <- 
  summary(fit_fs1)$fixed %>% 
  rownames_to_column(var= "term") %>% 
  filter(term!= "Intercept") %>% 
  mutate(predictor= "food sharing") %>% 
  mutate(sample= "phd bats")

# plot model results
(plot <- 
    coefs_fs1[1,] %>% 
    mutate(label= paste(predictor, "in", sample)) %>% 
    ggplot(aes(x=Estimate, y=label))+
    geom_vline(xintercept = 0, linetype= "solid", color= 'black')+
    geom_point(size=2)+
    geom_errorbarh(aes(xmin=`l-95% CI`, xmax=`u-95% CI`, height=0.1), size=1)+
    xlab("vocal distance")+
    ylab("")+
    theme_bw()+
    coord_cartesian(xlim= c(-1,1))+
    theme(axis.text=element_text(size=11), strip.text = element_text(size=12)))


# model dist as predicted by food sharing: Postdoc ------
fs_doc <- foodsharing %>% 
  filter(season1=="postdoc" & season2== "postdoc")

# MRQAP
fsm2 <- 
  fs_doc %>% 
  select(dyad, sharing) %>% 
  a_b_edgelist_to_matrix(directed=F )
fdm2 <- 
  fs_doc %>% 
  select(dyad, distance) %>% 
  a_b_edgelist_to_matrix(directed=F)
mrqap.dsp(fdm2~fsm2)  

# choose number of chains and chain length
# adjust as necessary
nchains <- 4
chain_length <- 3000
warmup_length <- 1000

# fit brms model of effect of familiarity on vocal distance
fit_fs2 <-
  brm(scale(distance) ~ 
        scale(sharing) + 
        (1|mm(actor,receiver)),
      data = fs_doc, 
      family = "gaussian",
      cores = nchains,
      chains = nchains,
      iter = chain_length,
      warmup = warmup_length)

# check model
hist(residuals(fit_fs2))
pp_check(fit_fs2, ndraws=100)
summary(fit_fs2)

# get model results
coefs_fs2 <- 
  summary(fit_fs2)$fixed %>% 
  rownames_to_column(var= "term") %>% 
  filter(term!= "Intercept") %>% 
  mutate(predictor= "food sharing") %>% 
  mutate(sample= "postdoc bats")

# plot model results
(plot <- 
    coefs_fs2[1,] %>% 
    mutate(label= paste(predictor, "in", sample)) %>% 
    ggplot(aes(x=Estimate, y=label))+
    geom_vline(xintercept = 0, linetype= "solid", color= 'black')+
    geom_point(size=2)+
    geom_errorbarh(aes(xmin=`l-95% CI`, xmax=`u-95% CI`, height=0.1), size=1)+
    xlab("vocal distance")+
    ylab("")+
    theme_bw()+
    coord_cartesian(xlim= c(-1,1))+
    theme(axis.text=element_text(size=11), strip.text = element_text(size=12)))

