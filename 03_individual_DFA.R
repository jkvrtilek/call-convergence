# linear discriminant function analysis to classify calls to bat
# Gerry Carter, gcarter1640@gmail.com
# updated by Julia Vrtilek 29 March 2024

# clear workspace
rm(list=ls())

# set working directory
setwd("/Users/jkvrtilek/Desktop/OSU/PhD/GitHub/call-convergence")

# load packages
library(MASS)
library(tidyverse)

# get data
raw <- readRDS("vampire_call_measures_filtered_transformed.RDS")

# get sample sizes
d <- raw %>% 
  separate(sound.files, into=c('ds','date', 'bat', 'file', 'call'), sep="_", remove = FALSE) %>% 
  group_by(bat) %>% 
  mutate(sample.size= n()) %>%
  ungroup()

# look at sample sizes
sort(unique(d$sample.size))

# set minimum sample size
# good rule of thumb is 20 obs per variable
predictors <- ncol(raw) - 3
min.sample.size <- predictors * 20

# filter calls by min sample size
d2 <- d %>% 
  filter(sample.size > min.sample.size)

# classify calls to bat using dfa with cross-validation (leave one-out classification)
dfa <- lda(bat ~ 
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
           CV=T, 
           data=d2)

# get classification matrix
cm <- table(d2$bat, dfa$class)

# get overall correct classification rate (accuracy)
# this is the best accuracy estimate
correct.cases <- sum(diag(cm))
all.cases <- sum(cm)
accuracy <- correct.cases/all.cases
accuracy

# get mean and range of correct assignments/all assignments to each bat
range(diag(cm)/colSums(cm), na.rm=T)
mean(diag(cm)/colSums(cm), na.rm=T)
# get mean and range of correct assignments/all calls from each bat
range(diag(cm)/rowSums(cm), na.rm=T)
mean(diag(cm)/rowSums(cm), na.rm=T)

# save classification matrix
write.csv(cm, "dfa-cv-results-standard.csv")




# plot accuracy as a function of sampling effort
tibble(bat= colnames(cm),
       accuracy= diag(cm)/colSums(cm),### I used column sums rather than row sums
       n.cases= colSums(cm)) %>% ###
  ggplot(aes(x=log(n.cases), y=accuracy))+
  geom_point(size=2)+
  geom_smooth(method= "lm")+
  xlab("log10 number of cases")+
  ylab('correct classification rate')+
  ggtitle('accuracy by sampling effort per bat')




# classify calls to bat using a single dfa without cross validation
# this method can give elevated accuracy (but doesn't matter because we are using this to get cross-classification rates)
dfa2 <- lda(bat ~ 
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
              peakf,
            CV=F, 
            data=d2)

# get classification rates
predictions <- predict(dfa2)
d2$prediction <- predictions$class
cm2 <- table(d2$bat, d2$prediction)

# get overall correct classification rate
correct.cases <- sum(diag(cm2))
all.cases <- sum(cm2)
accuracy2 <- correct.cases/all.cases
accuracy2

# get vocal distance as Mahalanobis distance between group centroids
meanGroup <- dfa2$means
distance <- as.matrix(dist(meanGroup %*% dfa2$scaling))

# save vocal similarity
write.csv(distance, file= "vocal-distance-lda.csv", row.names=F)





# get proportion of variance explained by the discriminant functions
props.df <- dfa2$svd*100/sum(dfa2$svd) 
tibble(prop= props.df, x= 1:length(props.df)) %>% 
  ggplot(aes(x=x, y=prop))+
  geom_col()+
  xlab("discriminant function")+
  ylab("proportion of variance explained")

# and the coefficients of the linear discriminant functions (these tell you which variables were most important for identifying cases)
dfa2$scaling

# save all DFA loadings sorted by absolute value of DF1
loadings <- 
  dfa2$scaling %>% 
  as.data.frame() %>% 
  arrange(desc(abs(LD1)))

write.csv(loadings, "dfa-loadings.csv")




# look at correlation between accuracy rates per bat (cross-validated or not)
n <- rowSums(cm)
c1 <- diag(cm)/n
c2 <- diag(cm2)/n

tibble(c1,c2, n) %>% 
  ggplot(aes(x=c2, y=c1))+
  geom_point(size=2)+
  geom_smooth(method="lm")+
  xlab("training accuracy")+
  ylab("testing accuracy")
cor.test(c1,c2)












# OPTIONAL: compare observed accuracy to random classification rates
# this takes awhile

# choose number of random datasets to test
perms <- 100

# store results in vector 'exp' for expected rate (by chance)
exp <-rep(NA, perms) 

# repeat dfa with random data
for (i in 1:perms) {
  
  # shuffle bat ids
  random.data <- 
    d4 %>% 
    mutate(bat= sample(bat))
  
  # classify calls to bat using dfa without cross validation
  dfa.rand <- lda(bat ~ 
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
                    peakf,
                  CV= F, 
                  data=random.data)
  
  # get classification rates
  predictions.rand <- predict(dfa.rand)
  d4$prediction.rand <- predictions.rand$class
  cm3 <- table(d4$prediction.rand, d4$bat)
  
  # get overall correct classification rate
  correct.cases <- sum(diag(cm3))
  all.cases <- sum(cm3)
  exp[i] <- correct.cases/all.cases
  
  # print progress
  print(paste(i, "of", perms))
}

# function to plot permutation test results----
hist_perm <- function(exp=exp, obs=obs, perms=perms, label=''){
  exp.range <- round(quantile(exp, probs= c(0.025, 0.975), na.rm=T),3)
  ggplot()+
    geom_histogram(aes(x=exp), color="black",fill="light blue")+
    geom_vline(aes(xintercept=obs), color="red", size=1)+
    xlab("expected values from null model")+
    ggtitle(label, subtitle = paste('obs = ',round(obs,3), ', exp 95% CI= ', exp.range[1], ' to ', exp.range[2], ", Prob exp >= obs: p", ifelse(mean(exp>=obs)==0,paste("<",1/perms), paste("=",signif(mean(exp>=obs),digits=2))),", permutations=",perms, sep=""))
}

# plot permutation test results
hist_perm(exp=exp, obs= accuracy2, perms=100)
