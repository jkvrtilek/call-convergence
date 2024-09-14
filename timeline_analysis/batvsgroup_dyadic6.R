# linear discriminant function analysis to classify calls: bat vs group
# Julia Vrtilek
# 12 Nov 2023

# load packages
library(tidyverse)
library(MASS)
library(graph4lg)
library(vegan)
library(boot)

# set file paths and wd
setwd("/Users/jkvrtilek/Desktop/OSU/PhD/GitHub/call-convergence/timeline_analysis")

# create functions
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


#### select relevant bats
# read in data
raw <- read_csv("metadata.csv")

# determine which bats we're using
# get tole and las pavas bats from 2016/2017
d <- raw %>% 
  distinct() %>% 
  filter(capture.site == "las.pavas" | capture.site == "tole") %>% 
  separate(date, into = c("year", "month", "day"), sep = "-", remove = FALSE) %>% 
  mutate(date = as.Date(date)) %>% 
  mutate(year = as.integer(year)) %>% 
  filter(year < 2018)

# get call numbers for each bat
sample <- read.csv("sample_size.csv")

# NOTE: chose to remove all bats with less than 100 calls!
las.pavas <- d %>%
  filter(capture.site == "las.pavas") %>%
  subset(bat_ID %in% sample$bat)

tole <- d %>%
  filter(capture.site == "tole") %>%
  subset(bat_ID %in% sample$bat)

# divide into groups by time
# add column of group names
d2 <- rbind(las.pavas,tole) %>% 
  mutate(
    group = case_when(capture.site == "las.pavas" ~ "lp",
                      capture.site == "tole" & date < "2016-08-01" ~ "tole1",
                      capture.site == "tole" & date > "2016-08-01" ~ "tole2")
  ) %>% 
  mutate(date.bat = paste(date,bat_ID,sep = "_"))


#### load and tidy spectral and temporal measurements
raw2 <- read.csv("/Users/jkvrtilek/Desktop/OSU/PhD/GitHub/LFS/vampire_call_measures_filtered_transformed.csv")

# select only relevant rows of specan
specan <- raw2 %>% 
  separate(sound.files, into = c("ds","date", "bat", "WAV", "sel"), sep = "_", remove = TRUE) %>% 
  dplyr::select(date:bat,duration:segments) %>% 
  mutate(date.bat = paste(date,bat,sep = "_")) %>% 
  subset(date.bat %in% d2$date.bat)

# select only relevant columns of bat metadata
temp <- d2 %>% 
  dplyr::select(group, date.bat)

# add group ID to measurements
specan2 <- left_join(specan, temp, by = "date.bat")

# rearrange columns into required order
specan2 <- specan2 %>% 
  mutate(individual = bat) %>% 
  dplyr::select(group, individual, 3:(ncol(specan2)-2))


#### make dataframe for each group
# tole pre-intro bats
tole_pre <- specan2 %>% 
  filter(group == "tole1")
# tole post-intro bats
tole_post <- specan2 %>% 
  filter(group == "tole2")
# las pavas bats
lp <- specan2 %>% 
  filter(group == "lp")


#### get all PRE-INTRO dyadic distances
pre <- rbind(tole_pre,lp)

# classify calls to bat using a single dfa without cross validation
dfa_pre <- lda(individual ~ 
              duration+    
              meanfreq+    
              sd+          
              freq.median+ 
              freq.Q25+ 
              freq.Q75+    
              freq.IQR+    
              time.median+
              time.Q25+   
              time.Q75+    
              time.IQR+    
              skew+        
              kurt+
              sp.ent+      
              time.ent+  
              entropy+     
              sfm+         
              meandom+     
              mindom+     
              maxdom+      
              dfrange+    
              modindx+     
              startdom+    
              enddom+      
              dfslope+   
              meanpeakf+   
              peakf+
              maxslope+
              minslope+
              abs_minslope+
              pos_slopes+
              neg_slopes+
              turns+
              meanslope+
              segments,
            CV= F, 
            data=pre)

# get vocal distance as Mahalanobis distance between group centroids
pre_distance <- as.matrix(dist(dfa_pre$means %*% dfa_pre$scaling))

# save all DFA loadings sorted by absolute value of DF1
loadings <- 
  dfa_pre$scaling %>% 
  as.data.frame() %>% 
  arrange(desc(abs(LD1)))

write.csv(loadings, "pre_dfa-loadings.csv")



#### get all POST-INTRO dyadic distances
post <- rbind(tole_post,lp)

# classify calls to bat using a single dfa without cross validation
dfa_post <- lda(individual ~ 
              duration+    
              meanfreq+    
              sd+          
              freq.median+ 
              freq.Q25+ 
              freq.Q75+    
              freq.IQR+    
              time.median+
              time.Q25+   
              time.Q75+    
              time.IQR+    
              skew+        
              kurt+
              sp.ent+      
              time.ent+  
              entropy+     
              sfm+         
              meandom+     
              mindom+     
              maxdom+      
              dfrange+    
              modindx+     
              startdom+    
              enddom+      
              dfslope+   
              meanpeakf+   
              peakf+
              maxslope+
              minslope+
              abs_minslope+
              pos_slopes+
              neg_slopes+
              turns+
              meanslope+
              segments,
            CV= F, 
            data=post)

# get vocal distance as Mahalanobis distance between group centroids
post_distance <- as.matrix(dist(dfa_post$means %*% dfa_post$scaling))

# save all DFA loadings sorted by absolute value of DF1
loadings <- 
  dfa_post$scaling %>% 
  as.data.frame() %>% 
  arrange(desc(abs(LD1)))

write.csv(loadings, "post_dfa-loadings.csv")


#### reorder matrices
order <- c(unique(tole$bat_ID), unique(las.pavas$bat_ID))
pre_dist <- reorder_mat(pre_distance, order)
post_dist <- reorder_mat(post_distance, order)


#### make matrix containing change in distance (convergence)
diff <- pre_dist - post_dist


#### make familiar vs introduced matrix
lp_bats <- unique(las.pavas$bat_ID)
tole_bats <- unique(tole$bat_ID)

a <- matrix(nrow = 8, ncol = 8)

col1 <- c(rep(0,times=4),rep(1,times=4))
col2 <- c(rep(1,times=4),rep(0,times=4))
v <- c(rep(col1,times=4),rep(col2,times=4))
a[,] <- v

diag(a) <- NA

rownames(a) <- c(tole_bats,lp_bats)
colnames(a) <- c(tole_bats,lp_bats)


#### Mantel test
mantel(diff, a, na.rm=T, method= 'spearman')


#### Gerry's code below
# get within-dyad convergence
dd <-  
  tibble(bat1=rownames(diff)[row(diff)], 
         bat2=colnames(diff)[col(diff)], 
         diff=c(diff)) %>% 
  filter(bat1 != bat2) %>% 
  mutate(site1= ifelse(bat1 %in% lp_bats, "lp", "tole")) %>% 
  mutate(site2= ifelse(bat2 %in% lp_bats, "lp", "tole")) %>% 
  mutate(dyad= ifelse(bat1<bat2, paste(bat1,bat2, sep="_"), paste(bat2,bat1, sep="_"))) %>% 
  mutate(familiar= site1== site2) %>% 
  group_by(dyad, familiar) %>% 
  summarize(diff= mean(diff)) %>% 
  ungroup()

# positive diff means they're converging
# get mean convergence
dd %>% 
  group_by(familiar) %>% 
  summarize(diff= mean(diff))

# what percentage of unfamiliar dyads had increased similarity?
dd %>% 
  mutate(increase= diff>0) %>% 
  group_by(familiar) %>% 
  summarize(prop.increase= mean(increase))

# ditto with 95% CIs, should be similar to previous numbers
dd %>% 
  mutate(increase= diff>0) %>% 
  boot_ci2(y = .$increase, x = .$familiar)

