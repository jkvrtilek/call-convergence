# test correlation between kinship and calls
# Gerry Carter, gcarter1640@gmail.com
# adapted 10 April 2024 by Julia Vrtilek 

# clear workspace
rm(list=ls())

# set working directory
setwd("/Users/jkvrtilek/Desktop/OSU/PhD/GitHub/call-convergence")

# load packages
library(tidyverse)
library(vegan)
library(igraph)
library(cowplot)
library(patchwork)
library(boot)

###################################################################

# functions

# get mean and 95% CI of values x via bootstrapping
boot_ci <- function(x, perms=5000, bca=F) {
  get_mean <- function(x, d) {
    return(mean(x[d]))
  }
  x <- as.vector(na.omit(x))
  mean <- mean(x)
  if(bca){
    boot <- boot.ci(boot(data=x, 
                         statistic=get_mean, 
                         R=perms, 
                         parallel = "multicore", 
                         ncpus = 4), 
                    type="bca")
    low <- boot$bca[1,4]
    high <- boot$bca[1,5] 
  }else{
    boot <- boot.ci(boot(data=x, 
                         statistic=get_mean, 
                         R=perms, 
                         parallel = "multicore", 
                         ncpus = 4), 
                    type="perc")
    low <- boot$perc[1,4]
    high <- boot$perc[1,5] 
  }
  c(low=low,mean=mean,high=high, N=round(length(x)))
}

# get mean and 95% CI via bootstrapping of values y within grouping variable x
boot_ci2 <- function(d=d, y=d$y, x=d$x, perms=5000, bca=F){
  df <- data.frame(effect=unique(x))
  df$low <- NA
  df$mean <- NA
  df$high <- NA
  df$n.obs <- NA
  for (i in 1:nrow(df)) {
    ys <- y[which(x==df$effect[i])]
    if (length(ys)>1 & var(ys)>0 ){
      b <- boot_ci(y[which(x==df$effect[i])], perms=perms, bca=bca) 
      df$low[i] <- b[1]
      df$mean[i] <- b[2]
      df$high[i] <- b[3]
      df$n.obs[i] <- b[4]
    }else{
      df$low[i] <- min(ys)
      df$mean[i] <- mean(ys)
      df$high[i] <- max(ys)
      df$n.obs[i] <- length(ys)
    }
  }
  df
}

# function convert matrix to edgelist
matrix_to_df <- function(m1){
  data.frame(dyad = paste(rownames(m1)[col(m1)], colnames(m1)[row(m1)], sep="_"),
             value = c(t(m1)), stringsAsFactors = FALSE)
}

# subset and sort common nodes in two networks (matrices)
# i.e. get the common individuals from two different networks to compare them
common_matrices <- function(m1 = m1, m2 = m2){
  # create subset of m1 for nodes that are in m2, and vice versa
  subset.m1 <- subset(m1, rownames(m1) %in% rownames(m2), 
                      select=c(colnames(m1) %in% colnames (m2)))
  subset.m2 <- subset(m2, rownames(m2) %in% rownames(subset.m1), 
                      select=c(colnames(m2) %in% colnames (subset.m1)))
  # put them in same order
  subset.m1 <- subset.m1[order(rownames(subset.m1)),order(colnames(subset.m1))]
  subset.m2 <- subset.m2[order(rownames(subset.m2)),order(colnames(subset.m2))]
  list(subset.m1,subset.m2)  
}

###################################################################

# get kinship data
k <- readRDS("vocal_social_data.RDS")

# get kinship matrix
km <- 
  k %>% 
  dplyr::select(bat1,bat2, kinship) %>% 
  graph_from_data_frame(directed=T) %>% 
  get.adjacency(attr='kinship', sparse=FALSE)

# get distances between group centroids
dist.m <- read.csv("vocal-distance-lda.csv") %>% 
  as.matrix() 
rownames(dist.m) <- colnames(dist.m)

# correlate kinship with vocal DISTANCE using mantel test
km2 <- common_matrices(km, dist.m)[[1]]
dm2 <- common_matrices(km, dist.m)[[2]]
set.seed(123)
mantel(km2, dm2, na.rm=T, method= 'spearman')

# add other social dyadic variables
d <- 
  k %>% 
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
    season_site == "2019_CH" ~ "same wild &\nsame lab colony",
    season_site == "2019_LB" ~ "same wild &\nsame lab colony",
    season_site == "2017_LP" ~ "same wild &\nsame lab colony",
    season_site == "2019_TL" ~ "same wild &\nsame lab colony",
    season_site == "2017_TL" ~ "same wild &\nsame lab colony",
    season_site == "2015_US" ~ "same long-term lab colony",
    season_site == "2017_different sites" ~ "different wild &\nsame lab colony",
    season_site == "2019_different sites" ~ "different wild &\nsame lab colony",
    season_site == "different colonies_TL" ~ "same wild &\ndifferent lab colony",
    season_site == "different colonies_different sites" ~ "different wild &\ndifferent lab colony",
    TRUE ~ "other"))

###################################################################

# plot vocal distance by time together
means <- 
  d %>% 
  filter(site!="other") %>% 
  filter(!is.na(dist)) %>% 
  boot_ci2(y=.$dist, x=.$treatment) %>% 
  rename(treatment = effect)

# plot means
p1 <- means %>% 
  mutate(treatment = fct_reorder(treatment, mean, .desc =T)) %>%
  filter(treatment != "same long-term lab colony") %>% 
  ggplot(aes(x=treatment, y=mean))+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=low, ymax=high, width=.1), size=1)+
  ylab("call distance")+
  xlab("")+
  ggtitle("vocal distance by time together")+
  theme_bw()+
  ylim(1.7, 2.7)
#+ theme(axis.text.y = element_text(size = 11), axis.title.y = element_text(size = 15))
p1

# plot vocal distance by kinship category
means_kin <- 
  d %>% 
  filter(!is.na(kinship)) %>% 
  mutate(kin= case_when(
    kinship == 0 ~ "nonkin",
    kinship >= 0 & kinship < 0.2 ~ "distant\nkin",
    kinship >= 0.2 ~ "close\nkin")) %>% 
  filter(!is.na(dist)) %>% 
  boot_ci2(y=.$dist, x=.$kin) %>% 
  rename(kin = effect)

# plot means
p2 <- 
  means_kin %>% 
  mutate(kin = fct_reorder(kin, mean, .desc =T)) %>%
  ggplot(aes(x=kin, y=mean))+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=low, ymax=high, width=.2), size=1)+
  ylab("call distance")+
  xlab("")+
  ggtitle("vocal distance by kinship")+
  theme_bw()+
  ylim(1.7, 2.7)
#+ theme(axis.text.y = element_text(size = 16), axis.text.x = element_text(size = 20), axis.title.y = element_text(size = 22))
p2

p2+p1+ plot_layout(widths = c(1, 4))

###################################################################

# violin plots

# kin
kinmeans <- means_kin %>% 
  mutate(kin = fct_reorder(kin, mean, .desc =F))

kindata <- d %>% 
  filter(!is.na(kinship)) %>% 
  mutate(kin = case_when(
    kinship == 0 ~ "nonkin",
    kinship >= 0 & kinship < 0.2 ~ "distant\nkin",
    kinship >= 0.2 ~ "close\nkin")) %>% 
  mutate(kin = ordered(kin, levels = c("nonkin", "distant\nkin", "close\nkin"))) %>%
  filter(!is.na(dist))

ggplot()+
  geom_point(data = kinmeans, aes(x=kin, y=mean), size=2)+
  geom_errorbar(data = kinmeans, aes(x=kin, y=mean, ymin=low, ymax=high, width=.2), size=1)+
  ylab("call distance")+
  xlab("")+
  geom_violin(data = kindata, aes(x=kin, y=dist), alpha = 0.1) +
  theme_bw()+
  #ylim(0.48, 0.65)+
  theme(axis.text.y = element_text(size = 20), axis.text.x = element_text(size = 24), axis.title.y = element_text(size = 30))

# familiarity
groupmeans <- means %>% 
  mutate(treatment = fct_reorder(treatment, mean, .desc =F)) %>%
  filter(treatment != "same long-term lab colony")

groupdata <- d %>% 
  filter(treatment != "same long-term lab colony") %>% 
  filter(treatment != "other") %>% 
  filter(!is.na(dist)) 

ggplot()+
  geom_point(data = groupmeans, aes(x=treatment, y=mean), size=2)+
  geom_errorbar(data = groupmeans, aes(x=treatment, y=mean, ymin=low, ymax=high, width=.1), size=1)+
  ylab("call distance")+
  xlab("")+
  geom_violin(data = groupdata, aes(x=treatment, y=dist), alpha = 0.1) +
  theme_bw()+
  #ylim(0.48, 0.65) +
  theme(axis.text.y = element_text(size = 16), axis.text.x = element_text(size = 24), axis.title.y = element_text(size = 30))


