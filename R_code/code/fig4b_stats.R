## Stats for figure 4B - comparing slopes of pRF density vs. eccentricity lines

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
             pattern = "*_radialFits_linear_contraOnly_20_40_10mmcontrol.mat")

for (f in files) {
  mf <- paste0(path,f)
  
  data <- readMat(mf)
  
  if (grepl("lh",mf)) {
    lh_rsquared = data$Adj.Rsquared
    lh_int = data$intercept
    lh_slope = data$slope
  } else if (grepl("rh",mf)) {
    rh_rsquared = data$Adj.Rsquared
    rh_int = data$intercept
    rh_slope = data$slope
  }
}

subject <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28)

# Organize right hemi data ---------------------------------------------------------------

hemi <- rep(c("right"), times = dim(rh_slope)[2])

stream <- rep(c("aventral"), times = dim(rh_slope)[2])
ROI <- rep(c("IOG"), times = dim(rh_slope)[2])
IOG <- data.frame(subject,hemi,stream,ROI,rh_slope[1,],rh_rsquared[1,])
colnames(IOG) <- c("subject","hemi","stream","ROI","slope","rsq")
ROI <- rep(c("pFus"), times = dim(rh_slope)[2])
pFus <- data.frame(subject,hemi,stream,ROI,rh_slope[2,],rh_rsquared[2,])
colnames(pFus) <- c("subject","hemi","stream","ROI","slope","rsq")
ROI <- rep(c("mFus"), times = dim(rh_slope)[2])
mFus <- data.frame(subject,hemi,stream,ROI,rh_slope[3,],rh_rsquared[3,])
colnames(mFus) <- c("subject","hemi","stream","ROI","slope","rsq")

stream <- rep(c("blateral"), times = dim(rh_slope)[2])
ROI <- rep(c("pSTS"), times = dim(rh_slope)[2])
pSTS <- data.frame(subject,hemi,stream,ROI,rh_slope[4,],rh_rsquared[4,])
colnames(pSTS) <- c("subject","hemi","stream","ROI","slope","rsq")
ROI <- rep(c("mSTS"), times = dim(rh_slope)[2])
mSTS <- data.frame(subject,hemi,stream,ROI,rh_slope[5,],rh_rsquared[5,])
colnames(mSTS) <- c("subject","hemi","stream","ROI","slope","rsq")

stream <- rep(c("blateral"), times = dim(rh_slope)[2])
ROI <- rep(c("CoS"), times = dim(rh_slope)[2])
CoS <- data.frame(subject,hemi,stream,ROI,rh_slope[6,],rh_rsquared[6,])
colnames(CoS) <- c("subject","hemi","stream","ROI","slope","rsq")

rh_org <- bind_rows(IOG, pFus, mFus, pSTS, mSTS,CoS)
levels(rh_org$ROI) <- c("IOG", "pFus", "mFus", "pSTS", "mSTS","CoS")
positions <- c("IOG", "pFus", "mFus", "pSTS", "mSTS","CoS")


# Right hemi stats -----------------------------------------------

# Some summary stats 
cdata <- ddply(rh_org, c("ROI","stream"), summarise,
               N    = length(slope),
               mean = mean(slope, na.rm=TRUE),
               se   = sem(slope),
               mean_rsq = mean(rsq, na.rm=TRUE),
               se_rsq   = sem(rsq)
)
cdata

# Organize for stats
rh_face <- bind_rows(IOG, pFus, mFus, pSTS, mSTS)
levels(rh_face$ROI) <- c("IOG", "pFus", "mFus", "pSTS", "mSTS")
rh_face$ROI <- factor(rh_face$ROI)
rh_face$stream <- factor(rh_face$stream)

# Subset data for paired t test
stream_means <- group_by(rh_face, subject, stream) %>%
  summarise_each(funs(mean(., na.rm = TRUE)), slope)
aventral <- subset(stream_means,  stream == "aventral", slope, drop=TRUE)
blateral <- subset(stream_means,  stream == "blateral", slope, drop=TRUE)
# Plot paired data (data viz only)
pd <- paired(aventral,blateral)
plot(pd, type = "profile") + theme_bw()

## Stats
# paired t-test by streams
rh_t <- t.test(aventral, blateral, paired = TRUE)
rh_t
## look by ROI
mod.lme <- lmer(slope ~ ROI + (1|subject), data = rh_face)
summary(mod.lme)
anova(mod.lme,type=c("III")) 
mod.lsm <- lsmeans::lsmeans(mod.lme, ~ROI)
mod.lsm
pairs <- contrast(mod.lsm,method="tukey")
pairs
pval <- summary(pairs)$p.value #exact sig figs for table

# Organize for hemi comp
rh_comp <- bind_rows(IOG, pFus, mFus, pSTS, mSTS)
levels(rh_comp$ROI) <- c("IOG", "pFus", "mFus", "pSTS", "mSTS")
rh_comp$ROI <- factor(rh_comp$ROI)
rh_comp$stream <- factor(rh_comp$stream)


# Organize left hemi data ---------------------------------------------------------------
hemi <- rep(c("left"), times = dim(lh_slope)[2])

stream <- rep(c("aventral"), times = dim(lh_slope)[2])
ROI <- rep(c("IOG"), times = dim(lh_slope)[2])
IOG <- data.frame(subject,hemi,stream,ROI,lh_slope[1,],lh_rsquared[1,])
colnames(IOG) <- c("subject","hemi","stream","ROI","slope","rsq")
ROI <- rep(c("pFus"), times = dim(lh_slope)[2])
pFus <- data.frame(subject,hemi,stream,ROI,lh_slope[2,],lh_rsquared[2,])
colnames(pFus) <- c("subject","hemi","stream","ROI","slope","rsq")
ROI <- rep(c("mFus"), times = dim(lh_slope)[2])
mFus <- data.frame(subject,hemi,stream,ROI,lh_slope[3,],lh_rsquared[3,])
colnames(mFus) <- c("subject","hemi","stream","ROI","slope","rsq")

stream <- rep(c("blateral"), times = dim(lh_slope)[2])
ROI <- rep(c("pSTS"), times = dim(lh_slope)[2])
pSTS <- data.frame(subject,hemi,stream,ROI,lh_slope[4,],lh_rsquared[4,])
colnames(pSTS) <- c("subject","hemi","stream","ROI","slope","rsq")
ROI <- rep(c("mSTS"), times = dim(lh_slope)[2])
mSTS <- data.frame(subject,hemi,stream,ROI,lh_slope[5,],lh_rsquared[55,])
colnames(mSTS) <- c("subject","hemi","stream","ROI","slope","rsq")

stream <- rep(c("blateral"), times = dim(lh_slope)[2])
ROI <- rep(c("CoS"), times = dim(lh_slope)[2])
CoS <- data.frame(subject,hemi,stream,ROI,lh_slope[6,],lh_rsquared[6,])
colnames(CoS) <- c("subject","hemi","stream","ROI","slope","rsq")

lh_org <- bind_rows(IOG, pFus, mFus, pSTS, mSTS,CoS)
levels(lh_org$ROI) <- c("IOG", "pFus", "mFus", "pSTS", "mSTS","CoS")
positions <- c("IOG", "pFus", "mFus", "pSTS", "mSTS","CoS")

## Stats
# Some summary stats 
cdata <- ddply(lh_org, c("ROI","stream"), summarise,
               N    = length(slope),
               mean = mean(slope, na.rm=TRUE),
               se   = sem(slope),
               mean_rsq = mean(rsq, na.rm=TRUE),
               se_rsq   = sem(rsq)
)
cdata

# Left hemi stats -----------------------------------------------

# Organize for stats
lh_face <- bind_rows(IOG, pFus, mFus, pSTS, mSTS)
levels(lh_face$ROI) <- c("IOG", "pFus", "mFus", "pSTS", "mSTS")
lh_face$ROI <- factor(lh_face$ROI)
lh_face$stream <- factor(lh_face$stream)

# Subset data for paired t test
lh_stream_means <- group_by(lh_face, subject, stream) %>%
  summarise_each(funs(mean(., na.rm = TRUE)), slope)
lh_aventral <- subset(lh_stream_means,  stream == "aventral", slope, drop=TRUE)
lh_blateral <- subset(lh_stream_means,  stream == "blateral", slope, drop=TRUE)
# Plot paired data
pd <- paired(lh_aventral,lh_blateral)
plot(pd, type = "profile") + theme_bw()

## Stats
# paired t-test by streams
lh_t <- t.test(lh_aventral, lh_blateral, paired = TRUE)
lh_t
## look by ROI
mod.lme <- lmer(slope ~ ROI + (1|subject), data = lh_face)
summary(mod.lme)
anova(mod.lme,type=c("III")) 
mod.lsm <- lsmeans::lsmeans(mod.lme, ~ROI)
mod.lsm
pairs <- contrast(mod.lsm,method="tukey")
pairs
pval <- summary(pairs)$p.value #exact sig figs for table

# Both hemis, face ROIs only
full <- full_join(rh_face, lh_face)
full$hemi <- factor(full$hemi)
full$stream <- factor(full$stream)
mod.lme <- lmer(slope ~ hemi*stream + (1|(subject)), data = full)
summary(mod.lme)
anova(mod.lme,type=c("III")) 
anova(mod.lme,type=c("III")) 