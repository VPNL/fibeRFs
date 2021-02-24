
rm(list=ls())
library(ggthemes)
library(R.matlab)
library(tidyverse)
library(lmerTest)
library(effectsize)
library(lsmeans)

# Load data ---------------------------------------------------------------
path <- "../results/study1/pRFs/"

files <- dir(paste0(path), 
             pattern = "*ratios.mat")

for (f in files) {
  mf <- paste0(path,f)
  
  data <- readMat(mf)
  
  if (grepl("lh",mf)) {
    lh_toon = data$ratioToon
    lh_wfret = data$ratioWfret
    
  } else if (grepl("rh",mf)) {
    rh_toon = data$ratioToon
    rh_wfret = data$ratioWfret
  }
}

subject <- c(1,2,3,4,5)

##Right hemi
hemi <- c("right","right","right","right","right")
# Toon
V1 <- rh_toon[,1]
IOG <- rh_toon[,2]
pFus <- rh_toon[,3]
mFus <- rh_toon[,4]
pSTS <- rh_toon[,5]
mSTS <- rh_toon[,6]
exper <- c("toon", "toon", "toon", "toon", "toon")
t <- data.frame(subject,hemi,exper,V1,IOG,pFus,mFus,pSTS,mSTS)

# Wfret
V1 <- rh_wfret[,1]
IOG <- rh_wfret[,2]
pFus <- rh_wfret[,3]
mFus <- rh_wfret[,4]
pSTS <- rh_wfret[,5]
mSTS <- rh_wfret[,6]
exper <- c("wfret", "wfret", "wfret", "wfret", "wfret")
w <- data.frame(subject,hemi,exper,V1,IOG,pFus,mFus,pSTS,mSTS)

rh <- full_join(t,w)

##Left hemi
hemi <- c("left","left","left","left","left")
# Toon
V1 <- lh_toon[,1]
IOG <- lh_toon[,2]
pFus <- lh_toon[,3]
mFus <- lh_toon[,4]
pSTS <- lh_toon[,5]
mSTS <- lh_toon[,6]
exper <- c("toon", "toon", "toon", "toon", "toon")
t <- data.frame(subject,hemi,exper,V1,IOG,pFus,mFus,pSTS,mSTS)

# Wfret
V1 <- lh_wfret[,1]
IOG <- lh_wfret[,2]
pFus <- lh_wfret[,3]
mFus <- lh_wfret[,4]
pSTS <- lh_wfret[,5]
mSTS <- lh_wfret[,6]
exper <- c("wfret", "wfret", "wfret", "wfret", "wfret")
w <- data.frame(subject,hemi,exper,V1,IOG,pFus,mFus,pSTS,mSTS)

lh <- full_join(t,w)

full <- full_join(lh,rh)

# Stats -----------------------------------------------------------------------
subs.tidy <- gather(full, ROI, mean, 4:9)
subs.tidy$ROI <- factor(subs.tidy$ROI)
subs.tidy$exper <- factor(subs.tidy$exper)
subs.tidy$hemi <- factor(subs.tidy$hemi)
mod.lme <- lmer(mean ~ hemi*exper*ROI + (1|subject), data = subs.tidy)
summary(mod.lme)
anova(mod.lme,type=c("III")) 
mod.lsm <- lsmeans::lsmeans(mod.lme, ~ ROI*exper)
mod.lsm
contrast(mod.lsm,method="tukey",by="ROI")

##By hemi separately
rh.tidy <- gather(rh, ROI, mean, 4:9)
rh.tidy$ROI <- factor(rh.tidy$ROI)
rh.tidy$exper <- factor(rh.tidy$exper)
rh.tidy$hemi <- factor(rh.tidy$hemi)
mod.lme <- lmer(mean ~ exper*ROI + (1|subject), data = rh.tidy)
summary(mod.lme)
ao <- anova(mod.lme,type=c("III")) 
ao #print anova results
#calculate effect size from test statistics
F_to_eta2(
  f = ao$`F value`,
  df = ao$NumDF,
  df_error = ao$DenDF
  )
#post-hoc tests
mod.lsm <- lsmeans::lsmeans(mod.lme, ~ ROI*exper)
mod.lsm
pairs<-contrast(mod.lsm,method="tukey",by="ROI")
pairs
pval <- summary(pairs)$p.value #exact p vals for table
t_to_d(t = summary(pairs)$t.ratio,
       df_error = summary(pairs)$df[1],
       pooled=TRUE) #get effect sizes

lh.tidy <- gather(lh, ROI, mean, 4:9)
lh.tidy$ROI <- factor(lh.tidy$ROI)
lh.tidy$exper <- factor(lh.tidy$exper)
lh.tidy$hemi <- factor(lh.tidy$hemi)
mod.lme <- lmer(mean ~ exper*ROI + (1|subject), data = lh.tidy)
summary(mod.lme)
ao <- anova(mod.lme,type=c("III")) 
ao #print anova results
#calculate effect size from test statistics
F_to_eta2(
  f = ao$`F value`,
  df = ao$NumDF,
  df_error = ao$DenDF
)
#post-hoc tests
mod.lsm <- lsmeans::lsmeans(mod.lme, ~ ROI*exper)
mod.lsm
pairs<-contrast(mod.lsm,method="tukey",by="ROI")
pairs
pval <- summary(pairs)$p.value #exact p vals for table
t_to_d(t = summary(pairs)$t.ratio,
       df_error = summary(pairs)$df[1],
       pooled=TRUE) #get effect sizes