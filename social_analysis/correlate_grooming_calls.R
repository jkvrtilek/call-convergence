# predict vocal similarity from grooming
# adapted 10 April 2024 from food-sharing code

library(tidyverse)
library(vegan)
library(asnipe)
library(boot)

# clear workspace
rm(list=ls())

# set working directory
setwd("/Users/jkvrtilek/Desktop/OSU/PhD/GitHub/call-convergence")

##################################################################

# function to convert rates to matrix
# convert a list of dyadic interactions to a network ------
a_b_edgelist_to_matrix <- function(el=el, symbol="_", directed=T, make.NA.zero=T){
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

# function convert matrix to edgelist
matrix_to_df <- function(m1){
  data.frame(dyad = paste(rownames(m1)[col(m1)], colnames(m1)[row(m1)], sep="_"),
             value = c(t(m1)), stringsAsFactors = FALSE)
}

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

# function to get mean and 95% CI via bootstrapping of values y within grouping variable x
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

# function to get the common individuals from two different networks to compare them
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

##################################################################

# get bond strength (grooming rates)
bonds <- read.csv("vocal_social_data.csv")

# use that function to get bond strength matrix (bm)
bm <- bonds %>% 
  dplyr::select(dir.dyad,grooming) %>% 
  filter(!is.na(grooming)) %>% 
  a_b_edgelist_to_matrix()

# get distances between group centroids
dist.m <- read.csv("vocal-distance-lda.csv") %>% 
  as.matrix() 
rownames(dist.m) <- colnames(dist.m)

# get kinship matrix
km <- bonds %>% 
  dplyr::select(actor,receiver,kinship) %>% 
  graph_from_data_frame(directed=T) %>% 
  get.adjacency(attr='kinship', sparse=FALSE)

##################################################################

# get matrices that share nodes
bm2 <- common_matrices(bm, dist.m)[[1]]
dm2 <- common_matrices(bm, dist.m)[[2]]
km2 <- common_matrices(dm2, km)[[2]]

# get correlation between bond strength and call similarity
set.seed(123)
mantel(bm2, dm2, method= "spearman")

# predict call similarity using both social bond and kinship
set.seed(123)
mrqap.dsp(dm2 ~ scale(bm2) + scale(km2), directed = "directed", test.statistic = "beta") 

# bats used
bats <- rownames(bm2)
test <- bonds %>%
  filter(actor %in% bats) %>% 
  select(actor, site1, season1, a.sex) %>% 
  distinct()

##################################################################

# get nonkin dyads
nonkin <- 
  matrix_to_df(km2) %>% 
  filter(value < 0.05) %>% 
  pull(dyad)

# add to bond data 
bonds2 <- 
  bonds %>% 
  filter(dir.dyad %in% nonkin)

# get means
means <- 
  bonds2 %>% 
  filter(!is.na(grooming)) %>% 
  mutate(categories = case_when(
    grooming == 0 ~ "no grooming\nobserved",
    grooming >= 0 ~ "grooming\nobserved")) %>% 
  filter(!is.na(dist)) %>% 
  boot_ci2(y=.$dist, x=.$categories) %>% 
  rename(categories = effect)

# get raw data
points <- 
  bonds2 %>% 
  filter(!is.na(grooming)) %>% 
  mutate(categories = case_when(
    grooming == 0 ~ "no grooming\nobserved",
    grooming >= 0 ~ "grooming\nobserved")) %>% 
  filter(!is.na(dist)) 

# plot means
p <- 
  means %>% 
  mutate(categories = fct_reorder(categories, mean, .desc =T)) %>%
  ggplot(aes(x=categories, y=mean))+
  geom_jitter(data= points, aes(y= dist), size=1, height=0, width=0.1, alpha=0.25)+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=low, ymax=high, width=.1), size=1)+
  ylab("vocal similarity")+
  xlab("")+
  ggtitle("nonkin call similarity by grooming relationship")+
  theme_bw()
p

# bats used
test2 <- bonds2 %>% 
  select(actor,site1,season1,a.sex) %>% 
  distinct()
