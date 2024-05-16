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

# functions
# convert a list of dyadic interactions to a network ------
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
d <- read_csv("/Users/gerry/Dropbox/Dropbox/_working/_ACTIVE/__students/Julia Vrtilek/2024/vocal_social_data.csv")

# need to fix food-sharing rates (I can see bat2 is never a donor which is wrong)
# need to label actor and receiver for directed dyad
# need to add sexes of both bats
# need to filter by sex and kinship

# model dist as predicted by food sharing
fs <- 
  d %>% 
  filter(!is.na(dist)) %>% 
  filter(!is.na(foodsharing)) %>% 
  filter(season1=="phd" & season2== "phd") %>% 
  mutate(dyad= ifelse(bat1<bat2, paste(bat1, bat2, sep= "_"), paste(bat2,bat1, sep= "_"))) %>% 
  arrange(dyad) %>% 
  group_by(dyad) %>% 
  summarize(kinship= mean(kinship, na.rm=T), sharing = mean(foodsharing, na.rm=T), distance= mean(dist, na.rm=T)) %>% 
  separate(dyad, into= c("bat1", "bat2"), sep= "_", remove=F) 
  
# MRQAP
fsm <- 
  fs %>% 
  select(dyad, sharing) %>% 
  a_b_edgelist_to_matrix(directed=F )
dm <- 
  fs %>% 
  select(dyad, distance) %>% 
  a_b_edgelist_to_matrix(directed=F)
mrqap.dsp(dm~fsm)  
  
# choose number of chains and chain length
# adjust as necessary
nchains <- 4
chain_length <- 3000
warmup_length <- 1000

# fit model of effect of familiarity on contact time
fit1 <-
  brm(scale(distance) ~ 
        scale(sharing) + 
        scale(kinship) +
        (1|mm(bat1,bat2)),
      data = fs, 
      family = "gaussian",
      cores = nchains,
      chains = nchains,
      iter = chain_length,
      warmup = warmup_length)

hist(residuals(fit1))
pp_check(fit1, ndraws=100)

summary(fit1)

# get model results
coefs1 <- 
  summary(fit1)$fixed %>% 
  rownames_to_column(var= "term") %>% 
  filter(term!= "Intercept") %>% 
  mutate(predictor= "food sharing") %>% 
  mutate(sample= "all bats")

# plot model results
(plot <- 
    coefs %>% 
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

  


