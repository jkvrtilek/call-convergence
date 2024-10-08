# Code for permuted DFA to compare calls before and after introduction
# Julia Vrtilek
# August 2023
# using Roger Mundry's pDFA script
# CHANGE FROM 6: cut off all bats with less than 100 calls

#### notes
# calls from same individual are not independent observations; therefore need permuted DFA
# 2 versions of permuted DFA: crossed and nested
# use nested: calls within individuals within groups. control factor = individual, grouping factor = group
# warning: in the script, there's a weird place where you input all the variables and it saves as a string and rewrites the model from a string

# could also assign calls to individual at three time steps and see whether individuals are maintaining individual signatures while converging - use normal DFA for this

#### setup
# load packages
library(MASS)
library(gridExtra)
library(grid)
library(lattice)
library(lemon)
library(tidyverse)

# set file paths and wd
setwd("/Users/jkvrtilek/Desktop/OSU/PhD/GitHub/call-convergence/timeline_analysis")

#### data wrangling
# read in data
raw <- read.csv("metadata.csv")

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
callnum <- read.csv("/Users/jkvrtilek/Desktop/OSU/PhD/GitHub/LFS/vampire_call_measures_filtered_transformed.csv") %>% 
  dplyr::select(sound.files) %>% 
  separate(sound.files, into = c("ds","date","bat","WAV","selec"), sep = "_") %>% 
  subset(bat %in% d$bat_ID) %>% 
  separate(date, into = c("year","month","day"), sep = "-", remove = FALSE) %>% 
  filter(year < 2018) %>% 
  filter(year > 2015) %>% 
  mutate(intro = case_when(year == 2016 & month < 8 ~ "pre",
                       year == 2016 & month > 8 ~ "post",
                       year == 2017 ~ "post")) %>% 
  group_by(bat, intro) %>%
  summarize(n = n()) %>% 
  pivot_wider(names_from = intro, values_from = n) %>% 
  filter(post > 100) %>% 
  arrange(post)

sample_size <- unique(callnum$bat)


# NOTE: chose to remove all bats with less than 100 calls!
las.pavas <- d %>%
  filter(capture.site == "las.pavas") %>%
  subset(bat_ID %in% sample_size)

tole <- d %>%
  filter(capture.site == "tole") %>%
  subset(bat_ID %in% sample_size)


#### divide into groups by time
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

# select only relevant rows
specan <- raw2 %>% 
  separate(sound.files, into = c("ds", "date", "bat", "WAV", "sel"), sep = "_", remove = TRUE) %>% 
  dplyr::select(date:bat,duration:segments) %>% 
  mutate(date.bat = paste(date,bat,sep = "_")) %>% 
  subset(date.bat %in% d2$date.bat)

# select only relevant columns
temp <- d2 %>% 
  dplyr::select(group, date.bat)

# add group ID to measurements
specan2 <- left_join(specan, temp, by = "date.bat")

# rearrange into required order
specan2 <- specan2 %>% 
  mutate(individual = bat) %>% 
  dplyr::select(group, individual, 3:(ncol(specan2)-2))


# Mundry's program requires the following:
# one column for test factor (groups)
# one column for control factor (individuals)
# however many variable columns
# variables should be only numbers, group and individual can be named
# NO missing values
# will take smallest number of calls

#for nested design, data do not have to be balanced
#tests for difference between groups ('testfac')
#groups and subjects ('contrfac') have to numbered consecutively and with integers beginning with 1

# make into function
pDFA <- function(xdata, test_fac, contr_fac, variables, n.sel = 100, nperm = 1000) {
  if (is.factor(test_fac)==F) {test_fac=as.factor(test_fac)}
  model=paste("lda(test_fac~",variables,", prior=pr_prob, subset=sel.index, data=xdata)",sep="")
  f.table=as.data.frame(table(contr_fac))
  ncf.levels=nrow(f.table)
  ntf.levels=nrow(table(test_fac))
  pr_prob=rep(1/ntf.levels, ntf.levels)#define prior probabilities to be equal for either of two groups
  #get number of cases per subject (level of contrfac):
  #get number of calls to select per subject (subject)
  n.to.sel=min(f.table$Freq)
  #set number of random selections original classification rate should be based on:
  ur.c.val=0
  ur.sel=0
  number=(1:nrow(xdata))
  #get assignment of subjects to groups
  gr=c()
  subj.gr= table(contr_fac,test_fac)
  subject=rownames(subj.gr)
  for (i in 1:ncf.levels){
    for (k in 1:ntf.levels){
      if (subj.gr[i,k]>0) {gr=c(gr,colnames(subj.gr)[k])}
    }
  }
  
  for (k in 1:n.sel){
    #make random selection of same number of cases per subject
    sel=rep(NA,nrow(xdata))#create var. for the random selection to be indicated
    for (i in 1:ncf.levels){
      sel[contr_fac==subject[i]] =sample(c(rep(1, n.to.sel), rep(0,f.table[i,2]-n.to.sel)), f.table[i,2],replace=F)
    }
    sel.index= number[sel==1]
    #do a DFA and store results in 'res':
    res=eval(parse(text=model))
    #get predictions and store them in 'pred':
    pred=predict(res,xdata,prior=pr_prob)$class
    ur.sel=ur.sel+sum((test_fac==pred)[sel==1])
    ur.c.val= ur.c.val+ sum((test_fac==pred)[sel==0])
  }
  #save number of correctly classified calls in variable 'orig.res':
  ur.sel= ur.sel/ n.sel
  ur.c.val= ur.c.val/ n.sel
  
  #set P-value to 1 (since original data should be treated as 1 permutation):
  p.sel=1
  p.c.val=1
  all.corr=matrix(NA, nrow=nperm, ncol=2)
  all.corr[1,1]=ur.sel
  all.corr[1,2]=ur.c.val
  
  if (length(gr)==ncf.levels){
    for (k in 1:(nperm-1)){
      #randomize subjects' assignments to groups:
      r.gr=sample(gr,length(gr), replace=F)
      for (i in 1:length(subject)){
        test_fac[contr_fac==subject[i]]=r.gr[i]
      }
      
      #make random selection or same number of cases per subject
      sel=rep(NA,nrow(xdata))#create var. for the random selection to be indicated
      for (i in 1:ncf.levels){
        sel[contr_fac==subject[i]] =sample(c(rep(1, n.to.sel), rep(0,f.table[i,2]-n.to.sel)), f.table[i,2],replace=F)
      }
      sel.index= number[sel==1]
      #do a DFA and store results in 'res':
      res=eval(parse(text=model))
      #get predictions and store them in 'pred':
      pred=predict(res,xdata,prior=pr_prob)$class
      ran.sel= sum((test_fac==pred)[sel==1])
      ran.c.val= sum((test_fac==pred)[sel==0])
      if (ran.sel>=ur.sel){p.sel = p.sel + 1}
      if (ran.c.val>= ur.c.val){p.c.val= p.c.val + 1}
      all.corr[k+1,1]=ran.sel
      all.corr[k+1,2]=ran.c.val
    }
    what=c("N correctly assigned, original, selected", "P for selected", "N correctly assigned, original, cross-validated", "P for cross-validated", "N groups (levels of test factor)", "N subjects (levels of control factor)", "N cases total", "N cases selected per subject","N selected total", "N permutations","N random selections")
    value=c(ur.sel,p.sel/nperm,ur.c.val,p.c.val/nperm,ntf.levels,ncf.levels,nrow(xdata),n.to.sel,n.to.sel*ncf.levels,nperm,n.sel)
    result=data.frame(what,value)
  }else{
    result="at least one subject is member of two groups; no test done"
  }
  result
  write.table(result,file=paste("/Users/jkvrtilek/Desktop/OSU/PhD/GitHub/call-convergence/timeline_analysis/pDFA",levels(test_fac)[1],"_",levels(test_fac)[2]),sep="",row.names=F,col.names=T)
  all.corr[,1]#comprises the number correctly classified selected objects for all
  #permutations (with the first value being that for the original data)
  all.corr[,2]#same for the cross-validated objects
  hist(all.corr[,1])#shows the frequency distribution
}


#### make dataframe for each group
# tole pre-intro bats
tole_pre <- specan2 %>% 
  filter(group == "tole1")
# tole post-intro bats
tole_post <- specan2 %>% 
  filter(group == "tole2")
# las pavas bats
lp_post <- specan2 %>% 
  filter(group == "lp")


#### 1
#### pre-introduction period Tole bats who are present at post-intro period 1
#### vs post-intro period 1 Las Pavas

xdata <- rbind(tole_pre, lp_post)
test_fac <- as.factor(as.character(xdata$group))
contr_fac <- as.factor(as.character(xdata$individual))
variables <- paste(colnames(specan2)[3:ncol(specan2)], collapse = "+")

pDFA(xdata, test_fac, contr_fac, variables)

# basic DFA for clustering plot
dfa <- lda(individual ~ 
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
           data=xdata)

# adapted from https://www.r-bloggers.com/2014/01/computing-and-visualizing-lda-in-r/
predictions <- predict(dfa)
dataset1 = data.frame(bat = xdata[,"individual"], group = xdata[,"group"], lda = predictions$x)

p1 <- ggplot(data = dataset1, aes(x = lda.LD1, y = lda.LD2, color = group)) + 
  geom_point(size = 0.1, alpha = 0.3) +
  labs(x = "LD1", y = "LD2", color = "") +
  scale_color_manual(values = c("#0f2d59","#3eb7c7"),
                     labels = c("Las Pavas",
                                "Tole")) +
  theme_minimal() +
  #theme(legend.position = "none") +
  stat_ellipse()

p1


#### 2
#### post-intro period 1 Tole bats
#### vs post-intro period 1 Las Pavas

xdata <- rbind(tole_post, lp_post)
test_fac <- as.factor(as.character(xdata$group))
contr_fac <- as.factor(as.character(xdata$individual))
variables <- paste(colnames(specan2)[3:ncol(specan2)], collapse = "+")

pDFA(xdata, test_fac, contr_fac, variables)

# basic DFA for clustering plot
dfa <- lda(individual ~ 
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
           data=xdata)

# adapted from https://www.r-bloggers.com/2014/01/computing-and-visualizing-lda-in-r/
predictions <- predict(dfa)
dataset2 = data.frame(bat = xdata[,"individual"], group = xdata[,"group"], lda = predictions$x)

p2 <- ggplot(data = dataset2, aes(x = lda.LD1, y = lda.LD2, color = group)) + 
  geom_point(size = 0.1, alpha = 0.3) +
  labs(x = "LD1", y = "LD2", color = "") +
  scale_color_manual(values = c("#0f2d59","#3eb7c7")) +
  theme_minimal() +
  theme(legend.position = "none") +
  stat_ellipse()

p2


subtitle1 <- textGrob("Post-introduction Las Pavas bats vs. pre-introduction Tole bats.")
subtitle2 <- textGrob("Post-introduction Las Pavas bats vs. post-introduction Tole bats.")
lay <- rbind(c(1,1),
             c(1,1),
             c(1,1),
             c(1,1),
             c(2,2),
             c(3,3),
             c(3,3),
             c(3,3),
             c(3,3),
             c(4,4))
grid1 <- grid_arrange_shared_legend(p1,subtitle1,p2,subtitle2,ncol=2,nrow=5,
                           top = "Title",layout_matrix = lay,position="bottom")
ggsave(paste("cluster_fig_", Sys.Date(), ".jpg", sep = ""),
       path = "/Users/jkvrtilek/Desktop/OSU/PhD/GitHub/call-convergence/results",
       plot = grid1)
