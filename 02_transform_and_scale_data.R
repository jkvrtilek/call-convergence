# script to transform data based on AD test, then scale it for DFA
# Julia Vrtilek
# 24 March 2024

# set working directory
setwd("/Users/jkvrtilek/Desktop/OSU/PhD/GitHub/call-convergence")

# load packages
library(tidyverse)
library(MASS)

# get data
raw <- readRDS("vampire_call_measures_filtered.RDS")

# make df to modify
d <- raw

# how to transform variables
wsqrt <- c("time.median","time.Q75","time.IQR","dfrange","abs_minslope","pos_slopes","neg_slopes")
wlog <- c("duration","time.Q25","skew","kurt","sfm","mindom","enddom","meanpeakf","turns")
winv <- c("freq.Q25","meandom","modindx")
rsqrt <- c("freq.median")
rlog <- c("freq.Q75","freq.IQR","segments")
rinv <- c("sp.ent","time.ent","entropy")
none <- c("meanfreq","sd","maxdom","startdom","dfslope","peakf","maxslope","minslope","meanslope")

# do transformations
# note: k = max + 1

# square root
for (i in 1:length(wsqrt)) {
  d[wsqrt[i]] <- sqrt(raw[wsqrt[i]])
}

# log
for (i in 1:length(wlog)) {
  d[wlog[i]] <- log(raw[wlog[i]]+1)
}

# inverse
for (i in 1:length(winv)) {
  d[winv[i]] <- 1/raw[winv[i]]
}

# reflect and root
for (i in 1:length(rsqrt)) {
  d[rsqrt[i]] <- sqrt(max(raw[rsqrt[i]])+1 - raw[rsqrt[i]])
}

# reflect and log
for (i in 1:length(rlog)) {
  d[rlog[i]] <- log(max(raw[rlog[i]],na.rm=T)+1 - raw[rlog[i]])
}

# reflect and inverse
for (i in 1:length(rinv)) {
  d[rinv[i]] <- 1/(max(raw[rinv[i]])+1-raw[rinv[i]])
}

# make function to fix scale() function so it returns a useful vector not a stupid matrix
scale2 <- function(x){as.vector(scale(x, scale = FALSE))}

# scale all numeric variables
d2 <- d %>% mutate(across(.cols=duration:peakf, .fns = scale2))

# filter out calls without fundfreq measures
d3 <- d2 %>% 
  filter(indicator == TRUE)

saveRDS(d3, "vampire_call_measures_filtered_transformed.RDS")




# # where indicator is FALSE, turn all fund freq measurements to 0
# ff <- colnames(d2)[30:37]
# 
# for (i in 1:nrow(d2)) {
#   if (d2[i,"indicator"] == FALSE) {
#     for (j in 1:length(ff)) {
#       d2[i,ff[j]] <- 0
#     }
#     print(i)
#   }
# }
