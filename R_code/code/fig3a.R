## Creates the basic plots for figure 3A (median size of pRFs)
# ROI indices are: 1 = V1, 2 = IOG, 3 = pFus, 4 = mFus, 5 = pSTS, 6 = mSTS, 7 = CoS

rm(list=ls())
library(R.matlab)
library(tidyverse)
library(lmerTest)
library(effectsize)
library(lsmeans)
library(dplyr)
library(PairedData)

sem <- function(x) {sd(x, na.rm=TRUE) / sqrt(sum(!is.na((x))))}

# Load data ---------------------------------------------------------------
path <- "../results/study1/pRFs/"

eccen <- 40
files <- dir(paste0(path), 
             pattern = "*_pRF_median_size_and_eccentricity_upTo40_ve20_10mmcontrol.mat")

for (f in files) {
  mf <- paste0(path,f)
  
  data <- readMat(mf)
  
  if (grepl("lh",mf)) {
    lh_sig = data$fix.sig
  } else if (grepl("rh",mf)) {
    rh_sig = data$fix.sig
  }
}

subject <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28)

# Organize right hemisphere -----------------------------------------------

hemi <- rep(c("right"), times = dim(rh_sig)[2])

stream <- rep(c("aventral"), times = dim(rh_sig)[2])
ROI <- rep(c("IOG"), times = dim(rh_sig)[2])
IOG <- data.frame(subject,hemi,stream,ROI,rh_sig[2,])
colnames(IOG) <- c("subject","hemi","stream","ROI","sig")
ROI <- rep(c("pFus"), times = dim(rh_sig)[2])
pFus <- data.frame(subject,hemi,stream,ROI,rh_sig[3,])
colnames(pFus) <- c("subject","hemi","stream","ROI","sig")
ROI <- rep(c("mFus"), times = dim(rh_sig)[2])
mFus <- data.frame(subject,hemi,stream,ROI,rh_sig[4,])
colnames(mFus) <- c("subject","hemi","stream","ROI","sig")

stream <- rep(c("blateral"), times = dim(rh_sig)[2])
ROI <- rep(c("pSTS"), times = dim(rh_sig)[2])
pSTS <- data.frame(subject,hemi,stream,ROI,rh_sig[5,])
colnames(pSTS) <- c("subject","hemi","stream","ROI","sig")
ROI <- rep(c("mSTS"), times = dim(rh_sig)[2])
mSTS <- data.frame(subject,hemi,stream,ROI,rh_sig[6,])
colnames(mSTS) <- c("subject","hemi","stream","ROI","sig")

stream <- rep(c("blateral"), times = dim(rh_sig)[2])
ROI <- rep(c("CoS"), times = dim(rh_sig)[2])
CoS <- data.frame(subject,hemi,stream,ROI,rh_sig[7,])
colnames(CoS) <- c("subject","hemi","stream","ROI","sig")

rh_org <- bind_rows(IOG, pFus, mFus, pSTS, mSTS,CoS)
levels(rh_org$ROI) <- c("IOG", "pFus", "mFus", "pSTS", "mSTS","CoS")
positions <- c("IOG", "pFus", "mFus", "pSTS", "mSTS","CoS")

# Plot right hemisphere ---------------------------------------------------

p <- ggplot(rh_org, aes(x=ROI, y=sig, col=ROI)) + 
  scale_x_discrete(limits = positions) +
  geom_boxplot(outlier.shape=NA) + #avoid plotting outliers twice
  geom_jitter(position=position_jitter(width=.1, height=0)) +
  labs(x = "Region of interest", y = "pRF size (in degrees)") +
  theme_classic() +
  theme(legend.position="none") +
  ylim(0,40)+
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values=c("#d73027", "#f46d43", "#fdae61", "#74add1","#4575b4", "#009900"))
  #scale_color_manual(values=c("#009900","#d73027", "#fdae61",  "#4575b4", "#f46d43",  "#74add1"))
p
ggsave("figs/3a_rh_40_ve20_10mmcontrol.pdf", p, width=7, height=5)

# Stats for right hemisphere ----------------------------------------------

stable <- ddply(rh_org, c("ROI"), summarise,
                N    = length(sig),
                mean = mean(sig, na.rm=TRUE),
                se   = sem(sig)
)
stable

#Organize for stats
rh_face <- bind_rows(IOG, pFus, mFus, pSTS, mSTS)
levels(rh_face$ROI) <- c("IOG", "pFus", "mFus", "pSTS", "mSTS")
rh_face$ROI <- factor(rh_face$ROI)
rh_face$stream <- factor(rh_face$stream)

#Subset data for paired t test
stream_means <- group_by(rh_face, subject, stream) %>%
                summarise_each(funs(mean(., na.rm = TRUE)), sig)
aventral <- subset(stream_means,  stream == "aventral", sig, drop=TRUE)
blateral <- subset(stream_means,  stream == "blateral", sig, drop=TRUE)
#Plot paired data for viz
pd <- paired(aventral,blateral)
plot(pd, type = "profile") + theme_bw()

#Paired t-test by streams
rh_t <- t.test(aventral, blateral, paired = TRUE)
rh_t
t_to_d(t = rh_t$statistic,
       df_error = rh_t$parameter,
       pooled=TRUE) #get effect sizes
#driven by pSTS?
mod.lme <- lmer(sig ~ ROI + (1|subject), data = rh_face)
summary(mod.lme)
ao <- anova(mod.lme,type=c("III")) 
ao #print anova results
#calculate effect size from test statistics
F_to_eta2(
  f = ao$`F value`,
  df = ao$NumDF,
  df_error = ao$DenDF
)
mod.lsm <- lsmeans::lsmeans(mod.lme, ~ROI)
mod.lsm
pairs<-contrast(mod.lsm,method="tukey")
pairs
pval <- summary(pairs)$p.value #exact p vals for table

#Organize data for later hemi comparison
rh_comp <- bind_rows(IOG, pFus, mFus, pSTS, mSTS)
levels(rh_comp$ROI) <- c("IOG", "pFus", "mFus", "pSTS", "mSTS")
rh_comp$ROI <- factor(rh_comp$ROI)
rh_comp$stream <- factor(rh_comp$stream)

# Organize left hemisphere -----------------------------------------------

hemi <- rep(c("left"), times = dim(lh_sig)[2])

stream <- rep(c("aventral"), times = dim(lh_sig)[2])
ROI <- rep(c("IOG"), times = dim(lh_sig)[2])
IOG <- data.frame(subject,hemi,stream,ROI,lh_sig[2,])
colnames(IOG) <- c("subject","hemi","stream","ROI","sig")
ROI <- rep(c("pFus"), times = dim(lh_sig)[2])
pFus <- data.frame(subject,hemi,stream,ROI,lh_sig[3,])
colnames(pFus) <- c("subject","hemi","stream","ROI","sig")
ROI <- rep(c("mFus"), times = dim(lh_sig)[2])
mFus <- data.frame(subject,hemi,stream,ROI,lh_sig[4,])
colnames(mFus) <- c("subject","hemi","stream","ROI","sig")

stream <- rep(c("blateral"), times = dim(lh_sig)[2])
ROI <- rep(c("pSTS"), times = dim(lh_sig)[2])
pSTS <- data.frame(subject,hemi,stream,ROI,lh_sig[5,])
colnames(pSTS) <- c("subject","hemi","stream","ROI","sig")
ROI <- rep(c("mSTS"), times = dim(lh_sig)[2])
mSTS <- data.frame(subject,hemi,stream,ROI,lh_sig[6,])
colnames(mSTS) <- c("subject","hemi","stream","ROI","sig")

stream <- rep(c("blateral"), times = dim(lh_sig)[2])
ROI <- rep(c("CoS"), times = dim(lh_sig)[2])
CoS <- data.frame(subject,hemi,stream,ROI,lh_sig[7,])
colnames(CoS) <- c("subject","hemi","stream","ROI","sig")

lh_org <- bind_rows(IOG, pFus, mFus, pSTS, mSTS,CoS)
levels(lh_org$ROI) <- c("IOG", "pFus", "mFus", "pSTS", "mSTS","CoS")
positions <- c("IOG", "pFus", "mFus", "pSTS", "mSTS","CoS")

# Plot left hemisphere ----------------------------------------------------

p <- ggplot(lh_org, aes(x=ROI, y=sig, col=ROI)) + 
  scale_x_discrete(limits = positions) +
  geom_boxplot(outlier.shape=NA) + #avoid plotting outliers twice
  geom_jitter(position=position_jitter(width=.1, height=0)) +
  labs(x = "Region of interest", y = "pRF size (in degrees)") +
  theme_classic() +
  theme(legend.position="none") +
  ylim(0,50)+
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values=c("#d73027", "#f46d43", "#fdae61", "#74add1","#4575b4", "#009900"))
p
ggsave("figs/3a_lh_40_ve20_10mmcontrol.pdf", p, width=7, height=5)

# Stats for left hemisphere -----------------------------------------------

stable <- ddply(lh_org, c("ROI"), summarise,
                N    = length(sig),
                mean = mean(sig, na.rm=TRUE),
                se   = sem(sig)
)
stable

lh_face <- bind_rows(IOG, pFus, mFus, pSTS, mSTS)
levels(lh_face$ROI) <- c("IOG", "pFus", "mFus", "pSTS", "mSTS")
lh_face$ROI <- factor(lh_face$ROI)
lh_face$stream <- factor(lh_face$stream)

# Subset data for paired t test
lh_stream_means <- group_by(lh_face, subject, stream) %>%
  summarise_each(funs(mean(., na.rm = TRUE)), sig)
lh_aventral <- subset(lh_stream_means,  stream == "aventral", sig, drop=TRUE)
lh_blateral <- subset(lh_stream_means,  stream == "blateral", sig, drop=TRUE)
# Plot paired data for vix
pd <- paired(lh_aventral,lh_blateral)
plot(pd, type = "profile") + theme_bw()

# paired t-test by streams
lh_t <- t.test(lh_aventral, lh_blateral, paired = TRUE)
lh_t
t_to_d(t = lh_t$statistic,
       df_error = lh_t$parameter,
       pooled=TRUE) #get effect sizes
## again driven by pSTS
mod.lme <- lmer(sig ~ ROI + (1|subject), data = lh_face)
summary(mod.lme)
ao <- anova(mod.lme,type=c("III")) 
ao #print anova results
#calculate effect size from test statistics
F_to_eta2(
  f = ao$`F value`,
  df = ao$NumDF,
  df_error = ao$DenDF
)
mod.lsm <- lsmeans::lsmeans(mod.lme, ~ROI)
mod.lsm
pairs<-contrast(mod.lsm,method="tukey")
pairs
pval <- summary(pairs)$p.value #exact p vals for table

# Both hemis, face ROIs only, no left mSTS
full <- full_join(rh_comp, lh_face)
full$hemi <- factor(full$hemi)
full$stream <- factor(full$stream)
mod.lme <- lmer(sig ~ hemi*stream + (1|(subject)), data = full)
summary(mod.lme)
anova(mod.lme,type=c("III")) 
mod.lsm <- lsmeans::lsmeans(mod.lme, ~ hemi*stream)
mod.lsm
contrast(mod.lsm,method="tukey") 