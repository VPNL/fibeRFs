## Creates figure 7B (fibers endpoints per eccentricity band)

rm(list=ls())
library(R.matlab)
library(tidyverse)
library(plyr)
library(lmerTest)
library(effectsize)
library(lsmeans)


sem <- function(x) {sd(x, na.rm=TRUE) / sqrt(sum(!is.na((x))))}


# Load data ---------------------------------------------------------------
path <- "../results/study1/fibers/10mm/"

files <- dir(paste0(path), 
             pattern = "*_weighted_4_bands.mat")

for (f in files) {
  mf <- paste0(path,f)
  
  data <- readMat(mf)
  
  if (grepl("lh",mf)) {
    lh_eccen = data$wprops
  } else if (grepl("rh",mf)) {
    rh_eccen = data$wprops
  }
}

subject <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)

# Organize right hemi ---------------------------------------------------------------
hemi <- rep(c("right"), times = dim(rh_eccen)[3])
# this is ugly but ¯\_(ツ)_/¯
stream <- rep(c("aventral"), times = dim(rh_eccen)[3])
ROI <- rep(c("IOG"), times = dim(rh_eccen)[3])
IOG <- data.frame(subject,hemi,stream,ROI,t(rh_eccen[1, , ]))
colnames(IOG) <- c("subject","hemi","stream","ROI",1,2,3,4)
IOG_gather <- gather(IOG,"bands","proportion", 5:8)
ROI <- rep(c("pFus"), times = dim(rh_eccen)[3])
pFus <- data.frame(subject,hemi,stream,ROI,t(rh_eccen[2, , ]))
colnames(pFus) <- c("subject","hemi","stream","ROI",1,2,3,4)
pFus_gather <- gather(pFus,"bands","proportion", 5:8)
ROI <- rep(c("mFus"), times = dim(rh_eccen)[3])
mFus <- data.frame(subject,hemi,stream,ROI,t(rh_eccen[3, , ]))
colnames(mFus) <- c("subject","hemi","stream","ROI",1,2,3,4)
mFus_gather <- gather(mFus,"bands","proportion", 5:8)

stream <- rep(c("blateral"), times = dim(rh_eccen)[3])
ROI <- rep(c("pSTS"), times = dim(rh_eccen)[3])
pSTS <- data.frame(subject,hemi,stream,ROI,t(rh_eccen[4, , ]))
colnames(pSTS) <- c("subject","hemi","stream","ROI",1,2,3,4)
pSTS_gather <- gather(pSTS,"bands","proportion", 5:8)
ROI <- rep(c("mSTS"), times = dim(rh_eccen)[3])
mSTS <- data.frame(subject,hemi,stream,ROI,t(rh_eccen[5, , ]))
colnames(mSTS) <- c("subject","hemi","stream","ROI",1,2,3,4)
mSTS_gather <- gather(mSTS,"bands","proportion", 5:8)
ROI <- rep(c("CoS"), times = dim(rh_eccen)[3])

stream <- rep(c("aventral"), times = dim(rh_eccen)[3])
CoS <- data.frame(subject,hemi,stream,ROI,t(rh_eccen[6, , ]))
colnames(CoS) <- c("subject","hemi","stream","ROI",1,2,3,4)
CoS_gather <- gather(CoS,"bands","proportion", 5:8)

#organize into one var
rh_eccen_org <- bind_rows(IOG_gather, pFus_gather, mFus_gather, pSTS_gather, mSTS_gather,CoS_gather)
levels(rh_eccen_org$ROI) <- c("IOG", "pFus", "mFus", "pSTS", "mSTS","CoS")
positions <- c("IOG", "pFus", "mFus", "pSTS", "mSTS","CoS")
#add alpha params for each band for nicer plotting
rh_eccen_org$alpha <- as.factor(ifelse(rh_eccen_org$bands ==1, 1,
                                       ifelse(rh_eccen_org$bands ==2,0.75,
                                              ifelse(rh_eccen_org$bands ==3,0.5,0.25))))

# Plot right hemisphere ---------------------------------------------------
b <- ggplot(rh_eccen_org, aes(x=ROI, y=proportion,col=bands,fill=ROI,alpha=factor(alpha))) +
  scale_x_discrete(limits = positions) +
  scale_y_continuous(limits = c(0,1)) +
  geom_bar(stat = "summary", fun.y = "mean",position = "dodge",width=.75,size=.3) +
  geom_point(aes(x=ROI),position=position_dodge(width=.75), alpha = 0.2,size = .5) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge",width=.75,alpha=.75,size=.3) +
  theme_classic()+
  labs(y = "Proportion fiber endpoints", x = "ROI") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_alpha_manual(values = c("0.25"=0.25, "0.5"=0.5, "0.75"=0.75, "1"=1), guide='none')+
  scale_fill_manual(values=c("#d73027", "#f46d43", "#fdae61", "#74add1","#4575b4", "#009900"), guide = 'none')+
  scale_colour_manual(values=c("#000000", "#000000", "#000000", "#000000", "#000000", "#000000"), guide = 'none') 
b
ggsave("figs/7b_rh.pdf", b, width=5, height=5)

# Stats for right hemisphere ---------------------------------------------------

# Some summary stats 
stable <- ddply(rh_eccen_org, c("ROI", "bands"), summarise,
               N    = length(proportion),
               mean = mean(proportion, na.rm=TRUE),
               se   = sem(proportion)
)
stable


# for stats
rh_face <- bind_rows(IOG_gather, pFus_gather, mFus_gather, pSTS_gather,mSTS_gather)
levels(rh_face$ROI) <- c("IOG", "pFus", "mFus", "pSTS", "mSTS")
rh_face$ROI <- factor(rh_face$ROI)
rh_face$stream <- factor(rh_face$stream)
mod.lme <- lmer(proportion ~ stream*bands + (1|subject), data = rh_face)
summary(mod.lme)
ao <- anova(mod.lme,type=c("III")) 
ao #print anova results
#calculate effect size from test statistics
F_to_eta2(
  f = ao$`F value`,
  df = ao$NumDF,
  df_error = ao$DenDF
)

mod.lsm <- lsmeans::lsmeans(mod.lme, ~ stream*bands)
mod.lsm
pairs <- contrast(mod.lsm,method="tukey",by="bands")
pairs
pval <- summary(pairs)$p.value #exact p vals 
t_to_d(t = summary(pairs)$t.ratio,
       df_error = summary(pairs)$df[1],
       pooled=TRUE) #get effect sizes


# Organize left hemi ---------------------------------------------------------------
hemi <- rep(c("left"), times = dim(lh_eccen)[3])
# this is ugly but ¯\_(ツ)_/¯
stream <- rep(c("aventral"), times = dim(lh_eccen)[3])
ROI <- rep(c("IOG"), times = dim(lh_eccen)[3])
IOG <- data.frame(subject,hemi,stream,ROI,t(lh_eccen[1, , ]))
colnames(IOG) <- c("subject","hemi","stream","ROI",1,2,3,4)
IOG_gather <- gather(IOG,"bands","proportion", 5:8)
ROI <- rep(c("pFus"), times = dim(lh_eccen)[3])
pFus <- data.frame(subject,hemi,stream,ROI,t(lh_eccen[2, , ]))
colnames(pFus) <- c("subject","hemi","stream","ROI",1,2,3,4)
pFus_gather <- gather(pFus,"bands","proportion", 5:8)
ROI <- rep(c("mFus"), times = dim(lh_eccen)[3])
mFus <- data.frame(subject,hemi,stream,ROI,t(lh_eccen[3, , ]))
colnames(mFus) <- c("subject","hemi","stream","ROI",1,2,3,4)
mFus_gather <- gather(mFus,"bands","proportion", 5:8)

stream <- rep(c("blateral"), times = dim(lh_eccen)[3])
ROI <- rep(c("pSTS"), times = dim(lh_eccen)[3])
pSTS <- data.frame(subject,hemi,stream,ROI,t(lh_eccen[4, , ]))
colnames(pSTS) <- c("subject","hemi","stream","ROI",1,2,3,4)
pSTS_gather <- gather(pSTS,"bands","proportion", 5:8)
ROI <- rep(c("mSTS"), times = dim(lh_eccen)[3])
mSTS <- data.frame(subject,hemi,stream,ROI,t(lh_eccen[5, , ]))
colnames(mSTS) <- c("subject","hemi","stream","ROI",1,2,3,4)
mSTS_gather <- gather(mSTS,"bands","proportion", 5:8)

stream <- rep(c("aventral"), times = dim(lh_eccen)[3])
ROI <- rep(c("CoS"), times = dim(lh_eccen)[3])
CoS <- data.frame(subject,hemi,stream,ROI,t(lh_eccen[6, , ]))
colnames(CoS) <- c("subject","hemi","stream","ROI",1,2,3,4)
CoS_gather <- gather(CoS,"bands","proportion", 5:8)

#organize into one var
lh_eccen_org <- bind_rows(IOG_gather, pFus_gather, mFus_gather, pSTS_gather, mSTS_gather,CoS_gather)
levels(lh_eccen_org$ROI) <- c("IOG", "pFus", "mFus", "pSTS", "mSTS","CoS")
positions <- c("IOG", "pFus", "mFus", "pSTS", "mSTS","CoS")
#add alpha params for each band for nicer plotting
lh_eccen_org$alpha <- as.factor(ifelse(lh_eccen_org$bands ==1, 1,
                                       ifelse(lh_eccen_org$bands ==2,0.75,
                                              ifelse(lh_eccen_org$bands ==3,0.5,0.25))))

# Plot left hemi ---------------------------------------------------------------
b <- ggplot(lh_eccen_org, aes(x=ROI, y=proportion,col=bands,fill=ROI,alpha=factor(alpha))) +
  scale_x_discrete(limits = positions) +
  scale_y_continuous(limits = c(0,1)) +
  geom_bar(stat = "summary", fun.y = "mean",position = "dodge",width=.75,size=.3) +
  geom_point(aes(x=ROI),position=position_dodge(width=.75), alpha = 0.2,size = .5) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge",width=.75,alpha=.75,size=.3) +
  theme_classic()+  
  labs(y = "Proportion fiber endpoints", x = "ROI") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_alpha_manual(values = c("0.25"=0.25, "0.5"=0.5, "0.75"=0.75, "1"=1), guide='none')+
  scale_fill_manual(values=c("#d73027", "#f46d43", "#fdae61", "#74add1","#4575b4", "#009900"), guide = 'none')+
  scale_colour_manual(values=c("#000000", "#000000", "#000000", "#000000", "#000000", "#000000"), guide = 'none') 
b
ggsave("figs/7b_lh.pdf", b, width=5, height=5)

# Stats for left hemi ---------------------------------------------------------------
# Some summary stats 
stable <- ddply(lh_eccen_org, c("ROI", "bands"), summarise,
                N    = length(proportion),
                mean = mean(proportion, na.rm=TRUE),
                se   = sem(proportion)
)
stable

# for stats
lh_face <- bind_rows(IOG_gather, pFus_gather, mFus_gather, pSTS_gather)
levels(lh_face$ROI) <- c("IOG", "pFus", "mFus", "pSTS")
lh_face$ROI <- factor(lh_face$ROI)
lh_face$stream <- factor(lh_face$stream)
mod.lme <- lmer(proportion ~ stream*bands + (1|subject), data = lh_face)
summary(mod.lme)
ao <- anova(mod.lme,type=c("III")) 
ao #print anova results
#calculate effect size from test statistics
F_to_eta2(
  f = ao$`F value`,
  df = ao$NumDF,
  df_error = ao$DenDF
)

mod.lsm <- lsmeans::lsmeans(mod.lme, ~ stream*bands)
mod.lsm
pairs <- contrast(mod.lsm,method="tukey",by="bands")
pairs
pval <- summary(pairs)$p.value #exact p vals 
t_to_d(t = summary(pairs)$t.ratio,
       df_error = summary(pairs)$df[1],
       pooled=TRUE) #get effect sizes