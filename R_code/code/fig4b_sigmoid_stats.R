rm(list=ls())
library(R.matlab)
library(tidyverse)
library(plyr)
library(lmerTest)
library(lsmeans)

sem <- function(x) {sd(x, na.rm=TRUE) / sqrt(sum(!is.na((x))))}

# Load data ---------------------------------------------------------------
path <- "../results/study1/pRFs/"

files <- dir(paste0(path), 
             pattern = "*_radialFits_sigmoid_contraOnly_20_40_10mmcontrol.mat")

# for stats
for (f in files) {
  mf <- paste0(path,f)
  
  data <- readMat(mf)
  
  if (grepl("lh",mf)) {
    lh_rsquared = data$Sig.Adj.Rsquared
    lh_max = data$Sig.a
    lh_min = data$Sig.b
    lh_50 = data$Sig.c
    lh_steep = data$Sig.d
    
  } else if (grepl("rh",mf)) {
    rh_rsquared = data$Sig.Adj.Rsquared
    rh_max = data$Sig.a
    rh_min = data$Sig.b
    rh_50 = data$Sig.c
    rh_steep = data$Sig.d

  }
}

subject <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28)

# Organize right hemi data ---------------------------------------------------------------

hemi <- rep(c("right"), times = dim(rh_50)[2])

stream <- rep(c("aventral"), times = dim(rh_50)[2])
ROI <- rep(c("IOG"), times = dim(rh_50)[2])
IOG <- data.frame(subject,hemi,stream,ROI,rh_50[1,],rh_min[1,],rh_max[1,],rh_steep[1,],rh_rsquared[1,])
colnames(IOG) <- c("subject","hemi","stream","ROI","fifty","min","max","steep","rsq")
ROI <- rep(c("pFus"), times = dim(rh_50)[2])
pFus <- data.frame(subject,hemi,stream,ROI,rh_50[2,],rh_min[2,],rh_max[2,],rh_steep[2,],rh_rsquared[2,])
colnames(pFus) <- c("subject","hemi","stream","ROI","fifty","min","max","steep","rsq")
ROI <- rep(c("mFus"), times = dim(rh_50)[2])
mFus <- data.frame(subject,hemi,stream,ROI,rh_50[3,],rh_min[3,],rh_max[3,],rh_steep[3,],rh_rsquared[3,])
colnames(mFus) <- c("subject","hemi","stream","ROI","fifty","min","max","steep","rsq")

stream <- rep(c("blateral"), times = dim(rh_50)[2])
ROI <- rep(c("pSTS"), times = dim(rh_50)[2])
pSTS <- data.frame(subject,hemi,stream,ROI,rh_50[4,],rh_min[4,],rh_max[4,],rh_steep[4,],rh_rsquared[4,])
colnames(pSTS) <- c("subject","hemi","stream","ROI","fifty","min","max","steep","rsq")
ROI <- rep(c("mSTS"), times = dim(rh_50)[2])
mSTS <- data.frame(subject,hemi,stream,ROI,rh_50[5,],rh_min[5,],rh_max[5,],rh_steep[5,],rh_rsquared[5,])
colnames(mSTS) <- c("subject","hemi","stream","ROI","fifty","min","max","steep","rsq")

stream <- rep(c("blateral"), times = dim(rh_50)[2])
ROI <- rep(c("CoS"), times = dim(rh_50)[2])
CoS <- data.frame(subject,hemi,stream,ROI,rh_50[6,],rh_min[6,],rh_max[6,],rh_steep[6,],rh_rsquared[6,])
colnames(CoS) <- c("subject","hemi","stream","ROI","fifty","min","max","steep","rsq")

rh_org <- bind_rows(IOG, pFus, mFus, pSTS, mSTS,CoS)
levels(rh_org$ROI) <- c("IOG", "pFus", "mFus", "pSTS", "mSTS","CoS")
positions <- c("IOG", "pFus", "mFus", "pSTS", "mSTS","CoS")



# Right hemi stats -----------------------------------------------

# Some summary stats 
cdata <- ddply(rh_org, c("ROI","stream"), summarise,
               N    = length(fifty),
               mean_50 = mean(fifty, na.rm=TRUE),
               se_50   = sem(fifty),
               mean_min = mean(min, na.rm=TRUE),
               se_min   = sem(min),
               mean_max = mean(max, na.rm=TRUE),
               se_max   = sem(max),
               mean_rsq = mean(rsq, na.rm=TRUE),
               se_rsq   = sem(rsq)
)
cdata

rh_face <- bind_rows(IOG, pFus, mFus, pSTS, mSTS)
levels(rh_face$ROI) <- c("IOG", "pFus", "mFus", "pSTS", "mSTS")
rh_face$ROI <- factor(rh_face$ROI)
rh_face$stream <- factor(rh_face$stream)

stream_means <- group_by(rh_face, subject, stream) %>%
  summarise_each(funs(mean(., na.rm = TRUE)), fifty, min, max, steep)

# Right hemi stats for 50pt-----------------------------------------------

# Subset data for paired t test
aventral <- subset(stream_means,  stream == "aventral", fifty, drop=TRUE)
blateral <- subset(stream_means,  stream == "blateral", fifty, drop=TRUE)

## Stats
# paired t-test by streams
rh_t <- t.test(aventral, blateral, paired = TRUE)
rh_t
## look by ROI
mod.lme <- lmer(fifty ~ ROI + (1|subject), data = rh_face)
summary(mod.lme)
anova(mod.lme,type=c("III")) 
mod.lsm <- lsmeans::lsmeans(mod.lme, ~ROI)
mod.lsm
contrast(mod.lsm,method="tukey")

# Plot paired data
pd <- paired(aventral,blateral)
plot(pd, type = "profile") + theme_bw()

# Right hemi stats for min pt-----------------------------------------------
aventral_min <- subset(stream_means,  stream == "aventral", min, drop=TRUE)
blateral_min <- subset(stream_means,  stream == "blateral", min, drop=TRUE)

## Stats
# paired t-test by streams
rh_t <- t.test(aventral_min, blateral_min, paired = TRUE)
rh_t
## look by ROI
mod.lme <- lmer(min ~ ROI + (1|subject), data = rh_face)
summary(mod.lme)
anova(mod.lme,type=c("III")) 
mod.lsm <- lsmeans::lsmeans(mod.lme, ~ROI)
mod.lsm
contrast(mod.lsm,method="tukey")

# Plot paired data
pd <- paired(aventral_min,blateral_min)
plot(pd, type = "profile") + theme_bw()

# Organize left hemi data ---------------------------------------------------------------

hemi <- rep(c("left"), times = dim(lh_50)[2])

stream <- rep(c("aventral"), times = dim(lh_50)[2])
ROI <- rep(c("IOG"), times = dim(lh_50)[2])
IOG <- data.frame(subject,hemi,stream,ROI,lh_50[1,],lh_min[1,],lh_max[1,],lh_steep[1,],lh_rsquared[1,])
colnames(IOG) <- c("subject","hemi","stream","ROI","fifty","min","max","steep","rsq")
ROI <- rep(c("pFus"), times = dim(lh_50)[2])
pFus <- data.frame(subject,hemi,stream,ROI,lh_50[2,],lh_min[2,],lh_max[2,],lh_steep[2,],lh_rsquared[2,])
colnames(pFus) <- c("subject","hemi","stream","ROI","fifty","min","max","steep","rsq")
ROI <- rep(c("mFus"), times = dim(lh_50)[2])
mFus <- data.frame(subject,hemi,stream,ROI,lh_50[3,],lh_min[3,],lh_max[3,],lh_steep[3,],lh_rsquared[3,])
colnames(mFus) <- c("subject","hemi","stream","ROI","fifty","min","max","steep","rsq")

stream <- rep(c("blateral"), times = dim(lh_50)[2])
ROI <- rep(c("pSTS"), times = dim(lh_50)[2])
pSTS <- data.frame(subject,hemi,stream,ROI,lh_50[4,],lh_min[4,],lh_max[4,],lh_steep[4,],lh_rsquared[4,])
colnames(pSTS) <- c("subject","hemi","stream","ROI","fifty","min","max","steep","rsq")
ROI <- rep(c("mSTS"), times = dim(lh_50)[2])
mSTS <- data.frame(subject,hemi,stream,ROI,lh_50[5,],lh_min[5,],lh_max[5,],lh_steep[5,],lh_rsquared[5,])
colnames(mSTS) <- c("subject","hemi","stream","ROI","fifty","min","max","steep","rsq")

stream <- rep(c("blateral"), times = dim(lh_50)[2])
ROI <- rep(c("CoS"), times = dim(lh_50)[2])
CoS <- data.frame(subject,hemi,stream,ROI,lh_50[6,],lh_min[6,],lh_max[6,],lh_steep[6,],lh_rsquared[6,])
colnames(CoS) <- c("subject","hemi","stream","ROI","fifty","min","max","steep","rsq")

lh_org <- bind_rows(IOG, pFus, mFus, pSTS, mSTS,CoS)
levels(lh_org$ROI) <- c("IOG", "pFus", "mFus", "pSTS", "mSTS","CoS")
positions <- c("IOG", "pFus", "mFus", "pSTS", "mSTS","CoS")


# Left hemi stats -----------------------------------------------

# Some summary stats 
cdata <- ddply(lh_org, c("ROI","stream"), summarise,
               N    = length(fifty),
               mean_50 = mean(fifty, na.rm=TRUE),
               se_50   = sem(fifty),
               mean_min = mean(min, na.rm=TRUE),
               se_min   = sem(min),
               mean_max = mean(max, na.rm=TRUE),
               se_max   = sem(max),
               mean_rsq = mean(rsq, na.rm=TRUE),
               se_rsq   = sem(rsq)
)
cdata

lh_face <- bind_rows(IOG, pFus, mFus, pSTS, mSTS)
levels(lh_face$ROI) <- c("IOG", "pFus", "mFus", "pSTS", "mSTS")
lh_face$ROI <- factor(lh_face$ROI)
lh_face$stream <- factor(lh_face$stream)

stream_means <- group_by(lh_face, subject, stream) %>%
  summarise_each(funs(mean(., na.rm = TRUE)), fifty, min, max, steep)

# Left hemi stats for 50pt-----------------------------------------------

# Subset data for paired t test
aventral <- subset(stream_means,  stream == "aventral", fifty, drop=TRUE)
blateral <- subset(stream_means,  stream == "blateral", fifty, drop=TRUE)

## Stats
# paired t-test by streams
lh_t <- t.test(aventral, blateral, paired = TRUE)
lh_t
## look by ROI
mod.lme <- lmer(fifty ~ ROI + (1|subject), data = lh_face)
summary(mod.lme)
anova(mod.lme,type=c("III")) 
mod.lsm <- lsmeans::lsmeans(mod.lme, ~ROI)
mod.lsm
contrast(mod.lsm,method="tukey")

# Plot paired data
pd <- paired(aventral,blateral)
plot(pd, type = "profile") + theme_bw()

# Left hemi stats for min pt-----------------------------------------------
aventral_min <- subset(stream_means,  stream == "aventral", min, drop=TRUE)
blateral_min <- subset(stream_means,  stream == "blateral", min, drop=TRUE)

## Stats
# paired t-test by streams
lh_t <- t.test(aventral_min, blateral_min, paired = TRUE)
lh_t
## look by ROI
mod.lme <- lmer(min ~ ROI + (1|subject), data = lh_face)
summary(mod.lme)
anova(mod.lme,type=c("III")) 
mod.lsm <- lsmeans::lsmeans(mod.lme, ~ROI)
mod.lsm
contrast(mod.lsm,method="tukey")

# Plot paired data
pd <- paired(aventral_min,blateral_min)
plot(pd, type = "profile") + theme_bw()
