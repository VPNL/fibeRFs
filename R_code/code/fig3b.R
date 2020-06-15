## Creates the basic plots for figure 3B (median size of pRFs in central eccen bands)
# Note that .mat files are organized such that 1st col is ROI, 2nd col is subject, 
# 3rd col is band #
# ROI indices are: 1 = V1, 2 = IOG, 3 = pFus, 4 = mFus, 5 = pSTS, 6 = mSTS, 7 = CoS

rm(list=ls())
library(ggthemes)
library(R.matlab)
library(tidyverse)

#set results path
path <- "../results/study1/pRFs/"

files <- dir(paste0(path), 
             pattern = "*_sigBandsControl_40_ve10.mat")

for (f in files) {
  mf <- paste0(path,f)
  
  data <- readMat(mf)
  
  if (grepl("lh",mf)) {
    lh_sig = data$sigByBand
  } else if (grepl("rh",mf)) {
    rh_sig = data$sigByBand
  }
}

subject <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21)

### Right hemisphere
V1 <- rh_sig[1, ,1]
IOG <- rh_sig[2, ,1]
pFus <- rh_sig[3, ,1]
mFus <- rh_sig[4, ,1]
pSTS <- rh_sig[5, ,1]
mSTS <- rh_sig[6, ,1]
CoS <- rh_sig[7, ,1]
rh <- data.frame(IOG,pFus,mFus,pSTS,mSTS,CoS)
rh.tidy_plot <- gather(rh, ROI, mean, 1:6)

rh_subs <- data.frame(subject,IOG,pFus,mFus,pSTS,mSTS,CoS)

positions <- c("IOG", "pFus", "mFus", "pSTS", "mSTS","CoS")

p <- ggplot(rh.tidy_plot, aes(x=ROI, y=mean)) + 
  scale_x_discrete(limits = positions) +
  geom_bar(aes(fill=ROI),stat = "summary", fun.y = "mean",position = "dodge",width=.75,size=.1) +
  geom_point(aes(fill=ROI),position=position_dodge(width=.75), alpha = 0.25,size = .75) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge",width=.75,alpha=.75,size=.3) +
  labs(x = "Region of interest", y = "Size (in degrees)") +
  theme_classic() +
  theme(legend.position="none") +
  ylim(0,40)+
  scale_fill_manual(values=c("#009900","#d73027", "#fdae61",  "#4575b4", "#f46d43",  "#74add1"))
p
ggsave("figs/3b_rh_40_ve10_band1.pdf", p, width=4, height=5)

## Stats
library(lmerTest)
rh_subs.tidy <- gather(rh_subs, ROI, mean, 2:6) #no CoS
rh_subs.tidy$subject <- factor(rh_subs.tidy$subject)
rh_subs.tidy$ROI <- factor(rh_subs.tidy$ROI)
mod.lme <- lmer(mean ~ ROI + (1|subject), data = rh_subs.tidy)
summary(mod.lme)
anova(mod.lme,type=c("III")) #specify type=c(“III”)to correct for the unbalanced design
library(lsmeans)
mod.lsm <- lsmeans::lsmeans(mod.lme, ~ ROI)
mod.lsm
pairs(mod.lsm)

## band 2
IOG <- rh_sig[2, ,2]
pFus <- rh_sig[3, ,2]
mFus <- rh_sig[4, ,2]
pSTS <- rh_sig[5, ,2]
mSTS <- rh_sig[6, ,2]
CoS <- rh_sig[7, ,2]
rh <- data.frame(IOG,pFus,mFus,pSTS,mSTS,CoS)
rh.tidy_plot <- gather(rh, ROI, mean, 1:6)

rh_subs <- data.frame(subject,IOG,pFus,mFus,pSTS,mSTS,CoS)

p <- ggplot(rh.tidy_plot, aes(x=ROI, y=mean)) + 
  scale_x_discrete(limits = positions) +
  geom_bar(aes(fill=ROI),stat = "summary", fun.y = "mean",position = "dodge",width=.75,size=.1) +
  geom_point(aes(fill=ROI),position=position_dodge(width=.75), alpha = 0.25,size = .75) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge",width=.75,alpha=.75,size=.3) +
  labs(x = "Region of interest", y = "Size (in degrees)") +
  theme_classic() +
  theme(legend.position="none") +
  ylim(0,40)+
  scale_fill_manual(values=c("#009900","#d73027", "#fdae61",  "#4575b4", "#f46d43",  "#74add1"))
p
ggsave("figs/3b_rh_40_ve10_band2.pdf", p, width=4, height=5)

## Stats
library(lmerTest)
rh_subs.tidy <- gather(rh_subs, ROI, mean, 2:6) #no CoS
rh_subs.tidy$subject <- factor(rh_subs.tidy$subject)
rh_subs.tidy$ROI <- factor(rh_subs.tidy$ROI)
mod.lme <- lmer(mean ~ ROI + (1|subject), data = rh_subs.tidy)
summary(mod.lme)
anova(mod.lme,type=c("III")) #specify type=c(“III”)to correct for the unbalanced design
library(lsmeans)
mod.lsm <- lsmeans::lsmeans(mod.lme, ~ ROI)
mod.lsm
pairs(mod.lsm)

## band 3
IOG <- rh_sig[2, ,3]
pFus <- rh_sig[3, ,3]
mFus <- rh_sig[4, ,3]
pSTS <- rh_sig[5, ,3]
mSTS <- rh_sig[6, ,3]
CoS <- rh_sig[7, ,3]
rh <- data.frame(IOG,pFus,mFus,pSTS,mSTS,CoS)
rh.tidy_plot <- gather(rh, ROI, mean, 1:6)

rh_subs <- data.frame(subject,IOG,pFus,mFus,pSTS,mSTS,CoS)

p <- ggplot(rh.tidy_plot, aes(x=ROI, y=mean)) + 
  scale_x_discrete(limits = positions) +
  geom_bar(aes(fill=ROI),stat = "summary", fun.y = "mean",position = "dodge",width=.75,size=.1) +
  geom_point(aes(fill=ROI),position=position_dodge(width=.75), alpha = 0.25,size = .75) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge",width=.75,alpha=.75,size=.3) +
  labs(x = "Region of interest", y = "Size (in degrees)") +
  theme_classic() +
  theme(legend.position="none") +
  ylim(0,40)+
  scale_fill_manual(values=c("#009900","#d73027", "#fdae61",  "#4575b4", "#f46d43",  "#74add1"))
p
ggsave("figs/3b_rh_40_ve10_band3.pdf", p, width=4, height=5)
