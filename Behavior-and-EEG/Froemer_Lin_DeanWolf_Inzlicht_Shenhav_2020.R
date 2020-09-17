library(R.matlab) # moving data between matlab and R (EEG)
library(lme4) # mixed effects models
library(MASS) # contrasts for LMMs
library(ggplot2) # Graphs
library(psych) # for corr.test
library(Hmisc) # for sedit when plotting topos
library(effects) # predicted effects
library(sjPlot) # tables
library(Rmisc) # data summaries and multiplot

# SPECIFY YOUR WORKING DIRECTORY! # 
setwd("./")

didLmerConverge = function(lmerModel){
  relativeMaxGradient=signif(max(abs(with(lmerModel@optinfo$derivs,solve(Hessian,gradient)))),3)
  if (relativeMaxGradient < 0.001) {
    cat(sprintf("\tThe relative maximum gradient of %s is less than our 0.001 criterion.\n\tYou can safely ignore any warnings about a claimed convergence failure.\n\n", relativeMaxGradient))
  }
  else {
    cat(sprintf("The relative maximum gradient of %s exceeds our 0.001 criterion.\nThis looks like a real convergence failure; maybe try simplifying your model?\n\n", relativeMaxGradient))
  }
}

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#########################  load and clean all data #################################

input_file = 'Behavior/FXC_102Data_010920.csv'
a1a = read.csv(input_file) # study 1 behavior only

a1a <- a1a[a1a$SubID > 1030,] # subjects with smaller IDs did an earlier Version of the paradigm without fixed reward rate

# behavior Study 2
input_file = 'Behavior/FXCallSubDataTable.csv'
a1 = read.csv(input_file)

## check participants performance an whether anyone needs to be excluded (cutoff 60% accuracy under high efficacy)
subAcc <- ddply(a1a, .(SubID, EffLvl), summarise,
                mresponse    = mean(Acc, na.rm=TRUE))

subAcc$SubID[subAcc$mresponse<0.6 & subAcc$EffLvl==1] # this applies to nobody in Study 1


## data prep for Study 1
a1a$SubID <-  as.factor(a1a$SubID)
a1a$RewLvl <- as.factor(a1a$RewLvl)
contrasts(a1a$RewLvl) <- contr.sdif(2)
a1a$EffLvl <- as.factor(a1a$EffLvl)
contrasts(a1a$EffLvl) <- contr.sdif(2)
a1a$Congruency <- as.factor(a1a$Congruency)
a1a$Congruency <- revalue(a1a$Congruency , c("0"="Incongruent", "1"="Congruent", "2"="Neutral"))
a1a$Congruency<- ordered(a1a$Congruency, levels=c("Incongruent", "Neutral","Congruent"))
contrasts(a1a$Congruency) <- contr.sdif(3) 
a1a$cTrial <- scale(a1a$Trial, scale=TRUE, center=TRUE)
a1a$fAcc <- as.factor(a1a$Acc)
contrasts(a1a$fAcc) <- contr.sdif(2)
a1a$IsRewarded <- as.factor(a1a$IsRewarded)
contrasts(a1a$IsRewarded) <- contr.sdif(2)

### data prep for Study 2
a1$SubID <-  as.factor(a1$SubID)
a1$RewLvl <- as.factor(a1$RewLvl)
contrasts(a1$RewLvl) <- contr.sdif(2)
a1$EffLvl <- as.factor(a1$EffLvl)
contrasts(a1$EffLvl) <- contr.sdif(2)
a1$Congruency <- as.factor(a1$Congruency)
a1$Congruency <- revalue(a1$Congruency , c("0"="Incongruent", "1"="Congruent", "2"="Neutral"))
a1$Congruency<- ordered(a1$Congruency, levels=c("Incongruent", "Neutral","Congruent"))
a1$cTrial <- scale(a1$Trial, scale=TRUE, center=TRUE)
contrasts(a1$Congruency) <- contr.sdif(3) #my.simple
a1$fAcc <- as.factor(a1$Acc)
contrasts(a1$fAcc) <- contr.sdif(2)
a1$IsRewarded <- as.factor(a1$IsRewarded)
contrasts(a1$IsRewarded) <- contr.sdif(2)

# get subset without missing responses for ERN analyses
a2 <- a1[!is.nan(a1$Resp),]
a2$Acc <- as.factor(a2$Acc)
contrasts(a2$Acc) <- contr.sdif(2)


for (i in levels(a1$SubID) ) 
{print(i)
  a1$zAccRT[a1$SubID==i] <-  scale(a1$AccRT[a1$SubID==i], scale= TRUE, center=TRUE )
  a1$zRT[a1$SubID==i] <-  scale(a1$RT[a1$SubID==i], scale= TRUE, center=TRUE )
}


for (i in levels(a1a$SubID) ) 
{print(i)
  a1a$zAccRT[a1a$SubID==i] <-  scale(a1a$AccRT[a1a$SubID==i], scale= TRUE, center=TRUE )
  a1a$zRT[a1a$SubID==i] <-  scale(a1a$RT[a1a$SubID==i], scale= TRUE, center=TRUE )
}


### EEG

input_file = 'Export/P3b250550.mat' # where is your file?
P3b = readMat(input_file)

P3b = P3b$P3b250550 #--> specify variable as saved in Matlab
P3b = as.data.frame(P3b) # convert to data frame
# specify column names
colnames(P3b)[c( 1, 2, 3, 4, 5, 
                 6, 7, 8, 9, 10, 
                 11, 12, 13, 14, 15, 
                 16, 17, 18, 19, 20, 
                 21, 22, 23, 24, 25, 
                 26, 27, 28, 29, 30,
                 31, 32, 33, 34, 35)]=c('Fp1', 'Fpz', 'Fp2', 'F8', 'F4', 
                                        'Fz', 'F3', 'F7', 'FC3', 'FCz', 
                                        'FC4', 'C4', 'Cz', 'C3', 'FT7', 
                                        'T3', 'FT8', 'T4', 'CP4', 'CPz', 
                                        'CP3', 'P3', 'Pz', 'P4', 'TP8', 
                                        'T6', 'TP7', 'T5', 'O1', 'Oz', 
                                        'O2', 'M1', 'M2', 'IO2', 'SO2')

a1$P3b <- apply(cbind (P3b$Pz,P3b$P4,P3b$P3), 1, mean) # computes average activation in ROI

input_file = 'Export/CNV10001500.mat' # where is your file?
CNV = readMat(input_file)

CNV = CNV$CNV10001500 #--> specify variable as saved in Matlab
CNV = as.data.frame(CNV) # convert to data frame
# specify column names
colnames(CNV)[c( 1, 2, 3, 4, 5, 
                 6, 7, 8, 9, 10, 
                 11, 12, 13, 14, 15, 
                 16, 17, 18, 19, 20, 
                 21, 22, 23, 24, 25, 
                 26, 27, 28, 29, 30,
                 31, 32, 33, 34, 35)]=c('Fp1', 'Fpz', 'Fp2', 'F8', 'F4', 
                                        'Fz', 'F3', 'F7', 'FC3', 'FCz', 
                                        'FC4', 'C4', 'Cz', 'C3', 'FT7', 
                                        'T3', 'FT8', 'T4', 'CP4', 'CPz', 
                                        'CP3', 'P3', 'Pz', 'P4', 'TP8', 
                                        'T6', 'TP7', 'T5', 'O1', 'Oz', 
                                        'O2', 'M1', 'M2', 'IO2', 'SO2')

a1$CNV <- apply(cbind (CNV$Fz,CNV$FCz,CNV$Cz), 1, mean) # computes average activation in ROI
# 
# input_file = '/Export/PreStimBl13001500.mat' # where is your file?
# PreSBl = readMat(input_file)
# 
# PreSBl = PreSBl$PreStimBl13001500 #--> specify variable as saved in Matlab
# PreSBl = as.data.frame(PreSBl) # convert to data frame
# # specify column names
# colnames(PreSBl)[c( 1, 2, 3, 4, 5, 
#                 6, 7, 8, 9, 10, 
#                 11, 12, 13, 14, 15, 
#                 16, 17, 18, 19, 20, 
#                 21, 22, 23, 24, 25, 
#                 26, 27, 28, 29, 30,
#                 31, 32, 33, 34, 35)]=c('Fp1', 'Fpz', 'Fp2', 'F8', 'F4', 
#                                        'Fz', 'F3', 'F7', 'FC3', 'FCz', 
#                                        'FC4', 'C4', 'Cz', 'C3', 'FT7', 
#                                        'T3', 'FT8', 'T4', 'CP4', 'CPz', 
#                                        'CP3', 'P3', 'Pz', 'P4', 'TP8', 
#                                        'T6', 'TP7', 'T5', 'O1', 'Oz', 
#                                        'O2', 'M1', 'M2', 'IO2', 'SO2')
# 
# a1$PreSBl <- apply(cbind (PreSBl$FC3, PreSBl$FCz, PreSBl$FC4, PreSBl$C3, PreSBl$C4, PreSBl$Cz), 1, mean) # computes average a


input_file = 'Export/Baseline2000.mat' # where is your file?
Baseline = readMat(input_file)

Baseline = Baseline$Baseline2000 #--> specify variable as saved in Matlab
Baseline = as.data.frame(Baseline) # convert to data frame
# specify column names
colnames(Baseline)[c( 1, 2, 3, 4, 5, 
                      6, 7, 8, 9, 10, 
                      11, 12, 13, 14, 15, 
                      16, 17, 18, 19, 20, 
                      21, 22, 23, 24, 25, 
                      26, 27, 28, 29, 30,
                      31, 32, 33, 34, 35)]=c('Fp1', 'Fpz', 'Fp2', 'F8', 'F4', 
                                             'Fz', 'F3', 'F7', 'FC3', 'FCz', 
                                             'FC4', 'C4', 'Cz', 'C3', 'FT7', 
                                             'T3', 'FT8', 'T4', 'CP4', 'CPz', 
                                             'CP3', 'P3', 'Pz', 'P4', 'TP8', 
                                             'T6', 'TP7', 'T5', 'O1', 'Oz', 
                                             'O2', 'M1', 'M2', 'IO2', 'SO2')

a1$Baseline <- apply(cbind (Baseline$CP4, Baseline$CPz, Baseline$CP3, Baseline$P3, Baseline$P4, Baseline$Pz), 1, mean) # computes average activation in ROI
a1$fBaseline <- apply(cbind (Baseline$Fz, Baseline$FCz, Baseline$Cz), 1, mean) # computes average activation in ROI


input_file = 'Export/rBaseline2000.mat' # where is your file?
rBaseline = readMat(input_file)

rBaseline = rBaseline$rBaseline2000 #--> specify variable as saved in Matlab
rBaseline = as.data.frame(rBaseline) # convert to data frame
# specify column names
colnames(rBaseline)[c( 1, 2, 3, 4, 5, 
                      6, 7, 8, 9, 10, 
                      11, 12, 13, 14, 15, 
                      16, 17, 18, 19, 20, 
                      21, 22, 23, 24, 25, 
                      26, 27, 28, 29, 30,
                      31, 32, 33, 34, 35)]=c('Fp1', 'Fpz', 'Fp2', 'F8', 'F4', 
                                             'Fz', 'F3', 'F7', 'FC3', 'FCz', 
                                             'FC4', 'C4', 'Cz', 'C3', 'FT7', 
                                             'T3', 'FT8', 'T4', 'CP4', 'CPz', 
                                             'CP3', 'P3', 'Pz', 'P4', 'TP8', 
                                             'T6', 'TP7', 'T5', 'O1', 'Oz', 
                                             'O2', 'M1', 'M2', 'IO2', 'SO2')

a2$rBaseline <- apply(cbind (rBaseline$FC4, rBaseline$FCz, rBaseline$FC3), 1, mean) # computes average activation in ROI

input_file = 'Export/P2PFCZ.mat' # where is your file?
P2PFRN = readMat(input_file)
a1$P2PFRN <- P2PFRN$P2PFCZ

input_file = 'Export/ERN0100.mat' # where is your file?
ERN = readMat(input_file)

ERN = ERN$ERN0100 #--> specify variable as saved in Matlab
ERN = as.data.frame(ERN) # convert to data frame
# specify column names
colnames(ERN)[c( 1, 2, 3, 4, 5, 
                 6, 7, 8, 9, 10, 
                 11, 12, 13, 14, 15, 
                 16, 17, 18, 19, 20, 
                 21, 22, 23, 24, 25, 
                 26, 27, 28, 29, 30,
                 31, 32, 33, 34, 35)]=c('Fp1', 'Fpz', 'Fp2', 'F8', 'F4', 
                                        'Fz', 'F3', 'F7', 'FC3', 'FCz', 
                                        'FC4', 'C4', 'Cz', 'C3', 'FT7', 
                                        'T3', 'FT8', 'T4', 'CP4', 'CPz', 
                                        'CP3', 'P3', 'Pz', 'P4', 'TP8', 
                                        'T6', 'TP7', 'T5', 'O1', 'Oz', 
                                        'O2', 'M1', 'M2', 'IO2', 'SO2')

a2$ERN <- apply(cbind (ERN$FC3, ERN$FCz, ERN$FC4), 1, mean) # computes average activation in ROI


# transformations
a1$P3bcor <- rep(NaN, length(a1$RT))
a1$CNVcor<- rep(NaN, length(a1$RT))
# cor --> residualized on baseline
a1$P3bcor[!is.nan(a1$P3b)] <- resid(lm(a1$P3b ~ a1$Baseline))
a1$CNVcor[!is.nan(a1$CNV)] <- resid(lm(a1$CNV ~ a1$fBaseline))

for (i in levels(a1$SubID) ) 
{print(i)
  
  a1$wsCNV[a1$SubID==i] <-  scale(a1$CNVcor[a1$SubID==i], scale= TRUE, center=TRUE )
  a1$wsP3b[a1$SubID==i] <-  scale(a1$P3bcor[a1$SubID==i], scale= TRUE, center=TRUE )
  a1$wsBaseline[a1$SubID==i] <-  scale(a1$Baseline[a1$SubID==i], scale= TRUE, center=TRUE )

}

## any subjects to exclude?
subAcc <- ddply(a1, .(SubID, EffLvl), summarise,
               mresponse    = mean(Acc, na.rm=TRUE))

subAcc$SubID[subAcc$mresponse<0.6 & subAcc$EffLvl==1]
a1o <- a1
## exclude thesubjects with accuracy in high efficacy lower than 60%
a1 <- a1[!a1$SubID==35 & !a1$SubID==1 & !a1$SubID==21 & !a1$SubID==38 & !a1$SubID==39 & !a1$SubID==43 & !a1$SubID==46 & !a1$SubID==53 & !a1$SubID==63,]
a1$SubID <-  as.factor(a1$SubID)

## also exclude subjects in EEG matrices for topo plotting later
P3bo <- P3b
CNVo <- CNV

P3b <- P3b[!a1o$SubID==35 & !a1o$SubID==1 & !a1o$SubID==21 & !a1o$SubID==38 & !a1o$SubID==39 & !a1o$SubID==43 & !a1o$SubID==46 & !a1o$SubID==53 & !a1o$SubID==63 , ]
CNV <- CNV[!a1o$SubID==35 & !a1o$SubID==1 & !a1o$SubID==21 & !a1o$SubID==38 & !a1o$SubID==39 & !a1o$SubID==43 & !a1o$SubID==46 & !a1o$SubID==53 & !a1o$SubID==63,]

## RT looks more or less normal
qplot(data=a1, x=RT, geom="density")+facet_grid(EffLvl~RewLvl) 

## code odd and even trials for reliability 
a1$oddeven <- a1$Trial %% 2
a1a$oddeven <- a1a$Trial %% 2

###### Success and Reward Rate Info ######
## table supplementary Method
## Reward Rate 
TimeoutsSubS1 <- ddply(a1a, .(SubID), summarise,
                       avpreDL    = length(IsRewarded[IsRewarded==1])/length(IsRewarded))

mean(TimeoutsSubS1$avpreDL)
sd(TimeoutsSubS1$avpreDL)

TimeoutsSubS2 <- ddply(a1, .(SubID), summarise,
                       avpreDL    = length(IsRewarded[IsRewarded==1])/length(IsRewarded))
mean(TimeoutsSubS2$avpreDL)
sd(TimeoutsSubS2$avpreDL)

## Success Rate 
TimeoutsSubS1 <- ddply(a1a, .(SubID), summarise,
                       avpreDL    = length(IsRewarded[Acc==1 & Fast==1])/length(IsRewarded))

mean(TimeoutsSubS1$avpreDL)
sd(TimeoutsSubS1$avpreDL)

TimeoutsSubS2 <- ddply(a1, .(SubID), summarise,
                       avpreDL    = length(IsRewarded[Acc==1 & Fast==1])/length(IsRewarded))

mean(TimeoutsSubS2$avpreDL)
sd(TimeoutsSubS2$avpreDL)


############################################################ Study 1 behavior only ############################################################
###############################################################################################################################################

## compute reliability for accuracy and RT
Reliability <- ddply(a1a[!a1a$IsMiss==1 & a1a$RT >200 ,], .(SubID, RewLvl, EffLvl, oddeven), summarise,
                     mAccRT    = mean(AccRT, na.rm=TRUE),
                     mAcc      = mean (Acc, na.rm=TRUE)
)

corr.test(Reliability$mAcc[Reliability$oddeven==1] ,Reliability$mAcc[Reliability$oddeven==0])
corr.test(Reliability$mAccRT[Reliability$oddeven==1] ,Reliability$mAccRT[Reliability$oddeven==0])

# Analyses:

# 1) How do reward and efficacy incentives affect behavior?

# 1a) Accuracy

# descriptives

sum1 <-summarySEwithin(a1a[!a1a$IsMiss,], measurevar="Acc", withinvars=c("EffLvl", "RewLvl"),  idvar="SubID", na.rm=TRUE)

pDRmCNVl3 <- ggplot(data=sum1, aes(x=EffLvl, y=Acc, group =RewLvl, fill=RewLvl))+geom_errorbar(aes(max = Acc+se, min = Acc-se), width=0.1,position=position_dodge(.1)) +  geom_point(shape = 22, size = 5, position=position_dodge(.1))+ scale_fill_manual(name = "Reward level", values=c("#343085" ,"#F1E923"))+
  theme_bw()+ geom_hline(aes(yintercept=0), colour="#000000", linetype=2, size=0.4)+ geom_line(position=position_dodge(.1))+ #facet_wrap(~ Congruency)+
  theme(axis.text = element_text(size = 12), strip.text.x= element_text(size=12))+xlab("Efficacy level") + ylab("Accuracy") + coord_cartesian(ylim=c(0.8,0.9)) +#ggtitle("Late CNV")+  #ylim(-0.12, 0.15)+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.position="none")

pdf("Figures/S1_Accuracynm_Efficacy_by_Reward.pdf", width =4, height = 4)
pDRmCNVl3
dev.off()

sum1 <-summarySEwithin(a1a[!a1a$IsMiss,], measurevar="Acc", withinvars=c("Congruency"),  idvar="SubID", na.rm=TRUE)

pDRmCNVl3 <- ggplot(data=sum1, aes(x=Congruency, y=Acc))+ geom_bar(stat="identity", position=position_dodge())+ #scale_fill_manual(name = "Reward level", values=c("#343085" ,"#F1E923"))+
  geom_errorbar(aes(max = Acc+se, min = Acc-se), width=0.1,position=position_dodge(.9)) + theme_bw()+ geom_hline(aes(yintercept=0), colour="#000000", linetype=2, size=0.4)+#facet_wrap(~ Congruency)+
  theme(axis.text = element_text(size = 12), strip.text.x= element_text(size=12))+xlab("Efficacy level") + ylab("Accuracy") + coord_cartesian(ylim=c(0.7,1)) +#ggtitle("Late CNV")+  #ylim(-0.12, 0.15)+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.position="none")

pdf("Figures/S1_Accuracynm_Congruency.pdf", width =4, height = 4)
pDRmCNVl3
dev.off()

print (summary(ACCmodnm <- glmer(Acc~ EffLvl*RewLvl + Congruency + cTrial+(Congruency|SubID), data = a1a[!a1a$IsMiss,], 
                                 family = binomial))) 

sv1_max <- svd(getME(ACCmodnm, "Tlist")[[1]]) 
sv1_max$d
round(sv1_max$d^2/sum(sv1_max$d^2)*100, 1)  

tab_model(ACCmodnm, transform = NULL)

didLmerConverge(ACCmodnm)

print (summary(ACCmodnmadd <- glmer(Acc~ EffLvl+RewLvl + Congruency + cTrial+(Congruency|SubID), data = a1a[!a1a$IsMiss,], 
                                 family = binomial))) 

sv1_max <- svd(getME(ACCmodnmadd, "Tlist")[[1]]) 
sv1_max$d
round(sv1_max$d^2/sum(sv1_max$d^2)*100, 1)  
didLmerConverge(ACCmodnmadd)
tab_model(ACCmodnmadd, transform = NULL)

anova(ACCmodnm, ACCmodnmadd)

# 1b) RT
# descriptives
sum1 <-summarySEwithin(a1a[!a1a$IsMiss & !a1a$RT<200,], measurevar="zAccRT", withinvars=c("EffLvl", "RewLvl"),  idvar="SubID", na.rm=TRUE)

pDRmCNVl <- ggplot(data=sum1, aes(x=EffLvl, y=zAccRT, group = RewLvl, fill=RewLvl))+geom_errorbar(aes(max = zAccRT+se, min = zAccRT-se), width=0.1,position=position_dodge(.1)) +  geom_point(shape = 22, size = 5, position=position_dodge(.1))+ scale_fill_manual(name = "Reward level", values=c("#343085" ,"#F1E923"))+
  theme_bw()+ geom_hline(aes(yintercept=0), colour="#000000", linetype=2, size=0.4)+ geom_line(position=position_dodge(.1))+ #facet_wrap(~ Congruency)+
  theme(axis.text = element_text(size = 12), strip.text.x= element_text(size=12))+xlab("Efficacy level") + ylab("Accurate RT (z-scored)") + coord_cartesian(ylim=c(-0.15,0.15))+ #ggtitle("Late CNV")+  #ylim(-0.12, 0.15)+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.position="none")

pdf("Figures/S1_AccurateZRT_Efficacy_by_Reward_dot_version.pdf", width =2, height = 4)
pDRmCNVl
dev.off()


sum1 <-summarySEwithin(a1a[!a1a$IsMiss & !a1a$RT<200,], measurevar="AccRT", withinvars=c("Congruency"),  idvar="SubID", na.rm=TRUE)

pDRmCNVl <- ggplot(data=sum1, aes(x=Congruency, y=AccRT))+ geom_bar(stat="identity", position=position_dodge())+ #scale_fill_manual(name = "Reward level", values=c("#343085" ,"#F1E923"))+
  geom_errorbar(aes(max = AccRT+se, min = AccRT-se), width=0.1,position=position_dodge(.9)) + theme_bw()+ geom_hline(aes(yintercept=0), colour="#000000", linetype=2, size=0.4)+ #facet_wrap(~ Congruency)+
  theme(axis.text = element_text(size = 12), strip.text.x= element_text(size=12))+xlab("Efficacy level") + ylab("Accurate RT (z-scored)") + coord_cartesian(ylim=c(550,700))+ #ggtitle("Late CNV")+  #ylim(-0.12, 0.15)+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.position="none")

pdf("Figures/S1_AccurateRT_Congruency.pdf", width =4, height = 4)
pDRmCNVl
dev.off()


print (summary(AccRTmod0t <- lmer(AccRT~ EffLvl*RewLvl +Congruency+cTrial+(Congruency|SubID), a1a[!a1a$IsMiss & !a1a$RT<200,], 
                                  REML=FALSE)))

sv1_max <- svd(getME(AccRTmod0t, "Tlist")[[1]]) 
sv1_max$d
round(sv1_max$d^2/sum(sv1_max$d^2)*100, 1)  

tab_model(AccRTmod0t)

qqnorm(resid(AccRTmod0t))
qqline(resid(AccRTmod0t))

tab_model(ACCmodnm,AccRTmod0t,transform=NULL) # table S1

###################################################### EEG Study ########################################################
#########################################################################################################################



Reliability <- ddply(a1[!a1$IsMiss==1 & a1$RT >200 ,], .(SubID, RewLvl, EffLvl, oddeven), summarise,
                     mAccRT    = mean(AccRT, na.rm=TRUE),
                     mAcc      = mean (Acc, na.rm=TRUE),
                     mP3b      = mean (P3bcor, na.rm=TRUE),
                     mCNV      = mean (CNVcor, na.rm=TRUE)
)



corr.test(Reliability$mAcc[Reliability$oddeven==1] ,Reliability$mAcc[Reliability$oddeven==0])
corr.test(Reliability$mAccRT[Reliability$oddeven==1] ,Reliability$mAccRT[Reliability$oddeven==0])
corr.test(Reliability$mP3b[Reliability$oddeven==1] ,Reliability$mP3b[Reliability$oddeven==0])
corr.test(Reliability$mCNV[Reliability$oddeven==1] ,Reliability$mCNV[Reliability$oddeven==0])

# Analyses:

# 1) How do reward and efficacy incentives affect behavior?

# 1a) Accuracy

# descriptives

sum1 <-summarySEwithin(a1[!a1$IsMiss,], measurevar="Acc", withinvars=c("EffLvl", "RewLvl"),  idvar="SubID", na.rm=TRUE)

pDRmCNVl3 <- ggplot(data=sum1, aes(x=EffLvl, y=Acc,group = RewLvl, fill=RewLvl))+geom_errorbar(aes(max = Acc+se, min = Acc-se), width=0.1,position=position_dodge(.1)) +  geom_point(shape = 22, size = 5, position=position_dodge(.1))+ scale_fill_manual(name = "Reward level", values=c("#343085" ,"#F1E923"))+
  theme_bw()+ geom_hline(aes(yintercept=0), colour="#000000", linetype=2, size=0.4)+ geom_line(position=position_dodge(.1))+ #facet_wrap(~ Congruency)+
  theme(axis.text = element_text(size = 12), strip.text.x= element_text(size=12))+xlab("Efficacy level") + ylab("Accuracy") + coord_cartesian(ylim=c(0.8,0.9)) +#ggtitle("Late CNV")+  #ylim(-0.12, 0.15)+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.position="none")

pdf("Figures/S2_Accuracynm_Efficacy_by_Reward.pdf", width =4, height = 4)
pDRmCNVl3
dev.off()

sum1 <-summarySEwithin(a1[!a1$IsMiss,], measurevar="Acc", withinvars=c("Congruency"),  idvar="SubID", na.rm=TRUE)

pDRmCNVl3 <- ggplot(data=sum1, aes(x=Congruency, y=Acc))+ geom_bar(stat="identity", position=position_dodge())+ #scale_fill_manual(name = "Reward level", values=c("#343085" ,"#F1E923"))+
  geom_errorbar(aes(max = Acc+se, min = Acc-se), width=0.1,position=position_dodge(.9)) + theme_bw()+ geom_hline(aes(yintercept=0), colour="#000000", linetype=2, size=0.4)+#facet_wrap(~ Congruency)+
  theme(axis.text = element_text(size = 12), strip.text.x= element_text(size=12))+xlab("Efficacy level") + ylab("Accuracy") + coord_cartesian(ylim=c(0.7,1)) +#ggtitle("Late CNV")+  #ylim(-0.12, 0.15)+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.position="none")

pdf("Figures/S2_Accuracynm_Congruency.pdf", width =4, height = 4)
pDRmCNVl3
dev.off()


print (summary(ACCmodnm <- glmer(Acc~ EffLvl*RewLvl + Congruency + cTrial+(Congruency|SubID), data = a1[!a1$IsMiss,], 
                                family = binomial))) 

sv1_max <- svd(getME(ACCmodnm, "Tlist")[[1]]) 
sv1_max$d
round(sv1_max$d^2/sum(sv1_max$d^2)*100, 1)  

tab_model(ACCmodnm, transform = NULL)

didLmerConverge(ACCmodnm)

# 1b) RT
# descriptives

sum1 <-summarySEwithin(a1[!a1$IsMiss & !a1$RT<200,], measurevar="zAccRT", withinvars=c("EffLvl", "RewLvl"),  idvar="SubID", na.rm=TRUE)

pDRmCNVl <- ggplot(data=sum1, aes(x=EffLvl, y=zAccRT, group = RewLvl, fill=RewLvl))+geom_errorbar(aes(max = zAccRT+se, min = zAccRT-se), width=0.1,position=position_dodge(.1)) +  geom_point(shape = 22, size = 5, position=position_dodge(.1))+ scale_fill_manual(name = "Reward level", values=c("#343085" ,"#F1E923"))+
  theme_bw()+ geom_hline(aes(yintercept=0), colour="#000000", linetype=2, size=0.4)+ geom_line(,position=position_dodge(.1))+ #facet_wrap(~ Congruency)+
  theme(axis.text = element_text(size = 12), strip.text.x= element_text(size=12))+xlab("Efficacy level") + ylab("Accurate RT (z-scored)") + coord_cartesian(ylim=c(-0.15,0.15))+ #ggtitle("Late CNV")+  #ylim(-0.12, 0.15)+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.position=c(0.5, 0.8))

pdf("Figures/S2_AccurateZRT_Efficacy_by_Reward_dot_version.pdf", width =2, height = 4)
pDRmCNVl
dev.off()


sum1 <-summarySEwithin(a1[!a1$IsMiss & !a1$RT<200,], measurevar="AccRT", withinvars=c("Congruency"),  idvar="SubID", na.rm=TRUE)

pDRmCNVl <- ggplot(data=sum1, aes(x=Congruency, y=AccRT))+ geom_bar(stat="identity", position=position_dodge())+ #scale_fill_manual(name = "Reward level", values=c("#343085" ,"#F1E923"))+
  geom_errorbar(aes(max = AccRT+se, min = AccRT-se), width=0.1,position=position_dodge(.9)) + theme_bw()+ geom_hline(aes(yintercept=0), colour="#000000", linetype=2, size=0.4)+ #facet_wrap(~ Congruency)+
  theme(axis.text = element_text(size = 12), strip.text.x= element_text(size=12))+xlab("Efficacy level") + ylab("Accurate RT (z-scored)") + coord_cartesian(ylim=c(550,700))+ #ggtitle("Late CNV")+  #ylim(-0.12, 0.15)+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.position="none")


pdf("Figures/S2_AccurateRT_Congruency.pdf", width =4, height = 4)
pDRmCNVl
dev.off()


print (summary(AccRTmod0t <- lmer(AccRT~ EffLvl*RewLvl +Congruency+cTrial+(EffLvl+Congruency|SubID), a1[!a1$IsMiss & !a1$RT<200,], 
                                  REML=FALSE)))

sv1_max <- svd(getME(AccRTmod0t, "Tlist")[[1]]) 
sv1_max$d
round(sv1_max$d^2/sum(sv1_max$d^2)*100, 1)  

tab_model(AccRTmod0t)

tab_model(ACCmodnm,AccRTmod0t,transform=NULL)

## as a control error trials only 

sum1 <-summarySEwithin(a1[!a1$IsMiss& a1$fAcc==0 & !a1$RT<200,], measurevar="zRT", withinvars=c("EffLvl", "RewLvl"),  idvar="SubID", na.rm=TRUE)

pDRmCNVl <- ggplot(data=sum1, aes(x=EffLvl, y=zRT, group = RewLvl, fill=RewLvl))+geom_errorbar(aes(max = zRT+se, min = zRT-se), width=0.1,position=position_dodge(.1)) +  geom_point(shape = 22, size = 5, position=position_dodge(.1))+ scale_fill_manual(name = "Reward level", values=c("#343085" ,"#F1E923"))+
  theme_bw()+ geom_hline(aes(yintercept=0), colour="#000000", linetype=2, size=0.4)+ geom_line(position=position_dodge(.1))+ #facet_wrap(~ Congruency)+
  theme(axis.text = element_text(size = 12), strip.text.x= element_text(size=12))+xlab("Efficacy level") + ylab("Error RT (z-scored)") + #ggtitle("Late CNV")+  #ylim(-0.12, 0.15)+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.position="none")


pdf("Figures/S2_ErrorRT_FXC_by_Reward.pdf", width =4, height = 4)
pDRmCNVl
dev.off()

print (summary(RTmod0Error <- lmer(RT~ (RewLvl*EffLvl+Congruency) +cTrial +(EffLvl+Congruency|SubID), a1[!a1$IsMiss & a1$fAcc==0 & !a1$RT<200,], 
                                       REML=FALSE)))

sv1_max <- svd(getME(RTmod0Error, "Tlist")[[1]]) 
sv1_max$d
round(sv1_max$d^2/sum(sv1_max$d^2)*100, 1)  

tab_model(RTmod0Error)
didLmerConverge(RTmod0Error)


# 2) How are incentive cues processed and how is reward allocated accordingly?
# 2a) P3b to cue
### P3b --> additive effects of efficacy and reward
print (summary(P3bmod0 <- lmer(P3bcor~ EffLvl*RewLvl+cTrial+Baseline+(EffLvl+RewLvl|SubID), a1[!a1$IsMiss & !a1$RT<200,], 
                               REML=FALSE)))

sv1_max <- svd(getME(P3bmod0, "Tlist")[[1]]) 
sv1_max$d
round(sv1_max$d^2/sum(sv1_max$d^2)*100, 1)  

tab_model(P3bmod0)

## Both Efficacy and Reward increase P3b amplitude
print (summary(P3bmod0t <- lmer(P3bcor ~ EffLvl+RewLvl+cTrial+Baseline+(EffLvl+RewLvl|SubID), a1[!a1$IsMiss & !a1$RT<200,], 
                                REML=FALSE)))

sv1_max <- svd(getME(P3bmod0t, "Tlist")[[1]]) 
sv1_max$d
round(sv1_max$d^2/sum(sv1_max$d^2)*100, 1)  

anova(P3bmod0, P3bmod0t)

sum1 <-summarySEwithin(a1[!a1$IsMiss & !a1$RT<200,], measurevar="P3b", withinvars=c("EffLvl", "RewLvl"),  idvar="SubID", na.rm=TRUE)

pDRmP3 <- ggplot(data=sum1, aes(x=EffLvl, y=P3b,group = RewLvl,  fill=RewLvl))+geom_errorbar(aes(max = P3b+se, min = P3b-se), width=0.1,position=position_dodge(.1)) +  geom_point(shape = 22, size = 5, position=position_dodge(.1))+ scale_fill_manual(name = "Reward level", values=c("#343085" ,"#F1E923"))+
  theme_bw()+  geom_line(position=position_dodge(.1))+ #facet_wrap(~ Congruency)+
  theme(axis.text = element_text(size = 12), strip.text.x= element_text(size=12))+xlab("Efficacy level") + ylab("Amplitude [µV]") + #ggtitle("Cue P3b")+  #ylim(-0.12, 0.15)+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.position="none")


pdf("Figures/P3_Efficacy_by_Reward.pdf", width =4, height = 4)
pDRmP3
dev.off()

# 2b) CNV

a1$cP3b <- scale(a1$P3bcor, scale=FALSE, center=TRUE) # P3b as covariate


sum1 <-summarySEwithin(a1[!a1$IsMiss & !a1$RT<200,], measurevar="CNV", withinvars=c("EffLvl", "RewLvl"),  idvar="SubID", na.rm=TRUE)

pDRmCNVl <- ggplot(data=sum1, aes(x=EffLvl, y=CNV,group = RewLvl, fill=RewLvl))+geom_errorbar(aes(max = CNV+se, min = CNV-se), width=0.1,position=position_dodge(.1)) +  geom_point(shape = 22, size = 5, position=position_dodge(.1))+ scale_fill_manual(name = "Reward level", values=c("#343085" ,"#F1E923"))+
  theme_bw()+  geom_line(position=position_dodge(.1))+ #facet_wrap(~ Congruency)+
  theme(axis.text = element_text(size = 12), strip.text.x= element_text(size=12))+xlab("Efficacy level") + ylab("Amplitude [µV]") + #ggtitle("Late CNV")+  #ylim(-0.12, 0.15)+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.position="none")

pdf("Figures/lCNV_Efficacy_by_Reward.pdf", width =4, height = 4)
pDRmCNVl
dev.off()


print (summary(CNVmod0 <- lmer(CNVcor~ (EffLvl*RewLvl)+(cP3b) +fBaseline+ cTrial+(EffLvl|SubID), a1[!a1$IsMiss & !a1$RT<200,], 
                               REML=FALSE)))

sv1_max <- svd(getME(CNVmod0, "Tlist")[[1]]) 
sv1_max$d
round(sv1_max$d^2/sum(sv1_max$d^2)*100, 1)  


print (summary(CNVmod0noERPs <- lmer(CNVcor~ (EffLvl*RewLvl) +fBaseline+ cTrial+(EffLvl|SubID), a1[!a1$IsMiss & !a1$RT<200,], 
                                     REML=FALSE)))

print (summary(CNVmod0noBL <- lmer(CNVcor~ (EffLvl*RewLvl)+(cP3b) + cTrial+(EffLvl|SubID), a1[!a1$IsMiss & !a1$RT<200,], 
                                     REML=FALSE)))

sv1_max <- svd(getME(CNVmod0noBL, "Tlist")[[1]]) 
sv1_max$d
round(sv1_max$d^2/sum(sv1_max$d^2)*100, 1)  


tab_model(P3bmod0, CNVmod0)

# 3) Do correlates of incentive processing predict performance on the upcoming trial?


# 3a) CNV and P3b on accuracy
# benchmark
print (summary(ACCmod0 <- glmer(Acc~ EffLvl*RewLvl + Congruency + cTrial+(Congruency|SubID), data = a1[!a1$IsMiss & !a1$RT<200,], 
                                family = binomial))) 

sv1_max <- svd(getME(ACCmod0, "Tlist")[[1]]) 
sv1_max$d
round(sv1_max$d^2/sum(sv1_max$d^2)*100, 1)  
## get residuals for control analysis 
a1$acmodRes <- rep(NA_real_, length(a1$SubID))
a1$acmodRes[!a1$IsMiss & !a1$RT<200] <- resid(ACCmod0)

## include ERPs
print (summary(ACCmod0P3bCNV <- glmer(Acc~ (EffLvl)*RewLvl+( wsP3b) +wsCNV+Congruency+ wsBaseline +cTrial+(Congruency|SubID), data = a1[!a1$IsMiss & !a1$RT<200,], 
                                      family = binomial))) 

sv1_max <- svd(getME(ACCmod0P3bCNV, "Tlist")[[1]]) 
sv1_max$d
round(sv1_max$d^2/sum(sv1_max$d^2)*100, 1)  

## compare effect sized of the two ERPs
library(car)
linearHypothesis(ACCmod0P3bCNV, c("-1*wsP3b - 1*wsCNV= 0"))

tab_model(ACCmod0P3bCNV, transform = NULL)

## plot predicted effect of ERPs on Accuracy
eff_df <- Effect(c( "wsCNV"), ACCmod0P3bCNV, xlevels=list(wsCNV =seq(min(a1$wsCNV, na.rm=TRUE ), max(a1$wsCNV, na.rm=TRUE), 0.1)))
#plot(eff_df)
IA <- as.data.frame(eff_df)
# here I revalue the levels of the factor so the legend is prettier
#IA$isChooseBestTrial<- revalue(IA$isChooseBestTrial, c("0"="worst", "1"="best"))
# now plot this with ggplot (this contains a bunch of nice things like shaded error bands, removing grids, changing axes labels etc.)
pAccRTCNV <- ggplot(data=IA, aes(x=wsCNV, y=fit )) +theme_bw(12) + geom_line()+
  geom_ribbon(data=IA, aes(x=wsCNV, max = fit + se, min = fit- se),alpha=0.1, inherit.aes = FALSE)+scale_y_continuous(limits=c(0.7,1), expand=c(0, 0))+ scale_x_continuous(limits = c(min(a1$wsCNV, na.rm=TRUE ), max(a1$wsP3b, na.rm=TRUE)),expand=c(0, 0))+
  xlab("CNV") + ylab("Accuracy") + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.position="bottom")

pdf("Figures/lCNV_acc.pdf", width =2, height = 4)
pAccRTCNV
dev.off()

eff_df <- Effect(c( "wsP3b"), ACCmod0P3bCNV, xlevels=list(wsP3b =seq(min(a1$wsP3b, na.rm=TRUE ), max(a1$wsP3b, na.rm=TRUE), 0.1)))
#plot(eff_df)
IA <- as.data.frame(eff_df)
# here I revalue the levels of the factor so the legend is prettier
#IA$isChooseBestTrial<- revalue(IA$isChooseBestTrial, c("0"="worst", "1"="best"))
# now plot this with ggplot (this contains a bunch of nice things like shaded error bands, removing grids, changing axes labels etc.)
pAccRTP3b <- ggplot(data=IA, aes(x=wsP3b, y=fit ))  +theme_bw(12) + geom_line()+
  geom_ribbon(data=IA, aes(x=wsP3b, max = fit + se, min = fit- se),alpha=0.1, inherit.aes = FALSE)+ scale_y_continuous(limits=c(0.7,1), expand=c(0, 0))+ scale_x_continuous(limits = c(min(a1$wsCNV, na.rm=TRUE ), max(a1$wsP3b, na.rm=TRUE)), expand=c(0, 0))+
  xlab("P3b") + ylab("Accuracy") + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.position="bottom")

pdf("Figures/P3b_acc.pdf", width =2, height = 4)
pAccRTP3b
dev.off()

##### test components only and get residuals for control analysis
print (summary(ACCmod0P3bCNVnoInc <- glmer(Acc~ ( wsP3b) +wsCNV+Congruency+ wsBaseline +cTrial+(Congruency|SubID), data = a1[!a1$IsMiss & !a1$RT<200,], 
                                      family = binomial))) 

sv1_max <- svd(getME(ACCmod0P3bCNVnoInc, "Tlist")[[1]]) 
sv1_max$d
round(sv1_max$d^2/sum(sv1_max$d^2)*100, 1)  

a1$acmodERPRes <- rep(NA_real_, length(a1$SubID))
a1$acmodERPRes[!a1$IsMiss & !a1$RT<200 & !is.nan(a1$wsP3b) & !is.nan(a1$wsCNV)] <- resid(ACCmod0P3bCNVnoInc)


### analyse residuals 
print (summary(ACCresmod0 <- lm(acmodERPRes~ EffLvl*RewLvl + Congruency + cTrial, data=a1[!a1$IsMiss & !a1$RT<200,], 
                                  REML=FALSE))) 

print (summary(ACCmod0ResP3bCNV <- lm(acmodRes~ (EffLvl)*RewLvl+( wsP3b) +wsCNV+Congruency+ wsBaseline +cTrial, data = a1[!a1$IsMiss & !a1$RT<200,], 
                                   REML=FALSE)))




# 3b) CNV and P3b on RT
## benchmark 

print (summary(AccRTmod0tb <- lmer(AccRT~ EffLvl*RewLvl +Congruency+cTrial+(Congruency|SubID), a1[!a1$IsMiss & !a1$RT<200,], 
                                   REML=FALSE)))

sv1_max <- svd(getME(AccRTmod0tb, "Tlist")[[1]]) 
sv1_max$d
round(sv1_max$d^2/sum(sv1_max$d^2)*100, 1)  

# get residuals for control analysis
a1$AccRTRes <- rep(NA_real_, length(a1$SubID))
a1$AccRTRes[!a1$IsMiss & !a1$RT<200 & !is.nan(a1$AccRT)] <- resid(AccRTmod0tb)

## with ERPs
print (summary(AccRTmod0t <- lmer(AccRT~ EffLvl*RewLvl +(wsP3b)+wsCNV +Congruency+ wsBaseline + cTrial+(EffLvl+wsCNV|SubID), a1[!a1$IsMiss & !a1$RT<200,],
                                  REML=FALSE)))

sv1_max <- svd(getME(AccRTmod0t, "Tlist")[[1]])
sv1_max$d
round(sv1_max$d^2/sum(sv1_max$d^2)*100, 1)

##
linearHypothesis(AccRTmod0t, c("-1*wsP3b - 1*wsCNV= 0"))

tab_model(ACCmod0P3bCNV, AccRTmod0t, transform = NULL)

## model w ERPs only - get residuals for control analysis
print (summary(AccRTmod0onlyERPs <- lmer(AccRT~  (wsP3b)+wsCNV +Congruency+ wsBaseline + cTrial+(wsCNV|SubID), a1[!a1$IsMiss & !a1$RT<200,],
                                         REML=FALSE)))

sv1_max <- svd(getME(AccRTmod0onlyERPs, "Tlist")[[1]])
sv1_max$d
round(sv1_max$d^2/sum(sv1_max$d^2)*100, 1)


a1$AccRTResNoInc <- rep(NA_real_, length(a1$SubID))

a1$AccRTResNoInc[!a1$IsMiss & !a1$RT<200 & !is.nan(a1$AccRT) & ! is.nan(a1$wsCNV) & ! is.nan(a1$wsP3b)] <- resid(AccRTmod0onlyERPs)

## analyze residuals

print (summary(AccRTREsmod0t <- lm(AccRTRes~ EffLvl*RewLvl +(wsP3b)+wsCNV +Congruency+ wsBaseline + cTrial, data =a1[!a1$IsMiss & !a1$RT<200,],
                                  REML=FALSE)))

tab_model(ACCmod0ResP3bCNV, AccRTREsmod0t)

## what do reward and efficacy do after we removed P3 and CNV
print (summary(AccRTREsmod0tnoERPs <- lm(AccRTResNoInc~ EffLvl*RewLvl +(wsP3b)+wsCNV +Congruency+ wsBaseline + cTrial, data=a1[!a1$IsMiss & !a1$RT<200,],
                                     REML=FALSE)))

sv1_max <- svd(getME(AccRTREsmod0tnoERPs, "Tlist")[[1]])
sv1_max$d
round(sv1_max$d^2/sum(sv1_max$d^2)*100, 1)

### plot predicted effects of ERPs on accurate RT

eff_df <- Effect(c( "wsCNV"), AccRTmod0t, xlevels=list(wsCNV =seq(min(a1$wsCNV, na.rm=TRUE ), max(a1$wsCNV, na.rm=TRUE), 0.1)))
#plot(eff_df)
IA <- as.data.frame(eff_df)
# here I revalue the levels of the factor so the legend is prettier
#IA$isChooseBestTrial<- revalue(IA$isChooseBestTrial, c("0"="worst", "1"="best"))
# now plot this with ggplot (this contains a bunch of nice things like shaded error bands, removing grids, changing axes labels etc.)
pAccRTCNV <- ggplot(data=IA, aes(x=wsCNV, y=fit )) +theme_bw(12) + geom_line()+
  geom_ribbon(data=IA, aes(x=wsCNV, max = fit + se, min = fit- se),alpha=0.1, inherit.aes = FALSE)+scale_y_continuous(limits=c(450,750), expand=c(0, 0))+ scale_x_continuous(limits = c(min(a1$wsCNV, na.rm=TRUE ), max(a1$wsP3b, na.rm=TRUE)),expand=c(0, 0))+
  xlab("CNV") + ylab("RT") + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.position="bottom")

pdf("Figures/lCNV_accRT.pdf", width =4, height = 4)
pAccRTCNV
dev.off()

eff_df <- Effect(c( "wsP3b"), AccRTmod0t, xlevels=list(wsP3b =seq(min(a1$wsP3b, na.rm=TRUE ), max(a1$wsP3b, na.rm=TRUE), 0.1)))
#plot(eff_df)
IA <- as.data.frame(eff_df)
# here I revalue the levels of the factor so the legend is prettier
#IA$isChooseBestTrial<- revalue(IA$isChooseBestTrial, c("0"="worst", "1"="best"))
# now plot this with ggplot (this contains a bunch of nice things like shaded error bands, removing grids, changing axes labels etc.)
pAccRTP3b <- ggplot(data=IA, aes(x=wsP3b, y=fit ))  +theme_bw(12) + geom_line()+
  geom_ribbon(data=IA, aes(x=wsP3b, max = fit + se, min = fit- se),alpha=0.1, inherit.aes = FALSE)+ scale_y_continuous(limits=c(450,750), expand=c(0, 0))+ scale_x_continuous(limits = c(min(a1$wsCNV, na.rm=TRUE ), max(a1$wsP3b, na.rm=TRUE)),expand=c(0, 0))+
  xlab("P3b") + ylab("RT") + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.position="bottom")

pdf("Figures/P3b_accRT.pdf", width =4, height = 4)
pAccRTP3b
dev.off()


# 5) How do incentives shape performance monitoring?

# 5a) ERN

a2 <- a2[!a2$SubID==35 & !a2$SubID==1 & !a2$SubID==21 & !a2$SubID==38 & !a2$SubID==39 & !a2$SubID==43 & !a2$SubID==46 & !a2$SubID==53 & !a2$SubID==63,]
a2$SubID <-  as.factor(a2$SubID)

a2$cRT <- a2$RT/1000
a2$crBaseline <- scale(a2$rBaseline, scale= FALSE, center= TRUE)

a2$numacc <- rep(0,length(a2$Acc))
a2$numacc[a2$Acc =="1"] <- 1


a2$meanErr <- rep(0,length(a2$Acc))

a2$TotErr <- rep(0,length(a2$Acc))
for (i in levels(a2$SubID) ) 
{print(i)
  
  for (j in levels(a2$EffLvl) ) {
    for (k in levels(a2$RewLvl) ) {
  a2$meanErr[a2$SubID==i & a2$EffLvl==j & a2$RewLvl==k] <-  mean(a2$numacc[a2$SubID==i & a2$EffLvl==j & a2$RewLvl==k] )
    }}
  a2$TotErr[a2$SubID==i] <-  mean(a2$numacc[a2$SubID==i] )
}

a2$meanCorr1 <- scale(a2$meanErr, scale=FALSE, center=TRUE)
a2$SubmeanCorr1 <- scale(a2$TotErr, scale=FALSE, center=TRUE)


print (summary(ERNmod0Con <- lmer(ERN~ (EffLvl*RewLvl*Acc)+Congruency*((cRT +I(cRT*cRT)))+(meanCorr1)*Acc +crBaseline+(Acc  +I(cRT*cRT)|SubID), a2[!a2$RT<200,], 
                                  REML=FALSE)))

sv1_max <- svd(getME(ERNmod0Con, "Tlist")[[1]]) 
sv1_max$d
round(sv1_max$d^2/sum(sv1_max$d^2)*100, 1)  


tab_model(ERNmod0Con)


didLmerConverge(ERNmod0Con)
anova(ERNmod0Con)

eff_df <- Effect(c("cRT", "Acc"), ERNmod0Con, xlevels=list(cRT =seq(min(a2$cRT, na.rm=TRUE ), max(a2$cRT, na.rm=TRUE), 0.01)))


IA <- as.data.frame(eff_df)
pRT<- ggplot(data=IA, aes(x=cRT, y=fit, color=Acc  )) + geom_line()+scale_colour_manual(name="Accuracy", values=cbPalette) +theme_bw(12)+ geom_ribbon(data=IA, aes(x=cRT, max = fit + se, min = fit- se, fill = Acc),alpha=0.1, inherit.aes = FALSE)+
  scale_fill_manual(name="Accuracy", values=cbPalette) + xlab("RT") + ylab("ERN") + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) #+ theme(legend.position=c(0.25, 0.2))

pdf("Figures/ERN_ACC_RT.pdf", width =4, height = 4)
pRT
dev.off()


eff_df <- Effect(c("meanCorr1", "Acc"), ERNmod0Con)
IA <- as.data.frame(eff_df)
pPI<- ggplot(data=IA, aes(x=meanCorr1, y=fit, color=Acc  )) + geom_line()+scale_colour_manual(name="Accuracy", values=cbPalette) +theme_bw(12)+ geom_ribbon(data=IA, aes(x=meanCorr1, max = fit + se, min = fit- se, fill = Acc),alpha=0.1, inherit.aes = FALSE)+
  scale_fill_manual(name="Accuracy", values=cbPalette) + xlab("P(Correct|Incentive)") + ylab("ERN") + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) #+ theme(legend.position=c(0.25, 0.2))

pdf("Figures/ERN_ACC_Exp.pdf", width =4, height = 4)
pPI
dev.off()

eff_df <- Effect(c("RewLvl", "EffLvl", "Acc"), ERNmod0Con)
IA <- as.data.frame(eff_df)
pRT<- ggplot(data=IA, aes(x=EffLvl, y=fit , group = interaction(RewLvl,Acc), fill= RewLvl, color = Acc)) +  geom_errorbar(aes(max = fit+se, min = fit-se), width=0.1,position=position_dodge(.1)) +geom_point(shape = 22, size = 5,position=position_dodge(.1))+geom_line()+ scale_fill_manual(name = "Reward level", values=c("#343085" ,"#F1E923"))+
   theme_bw()+ geom_hline(aes(yintercept=0), colour="#000000", linetype=2, size=0.4)+scale_colour_manual(name="Is Rewarded",values=c("#CC0000","#000000"))+#facet_wrap(~EffLvl)+
  theme(axis.text = element_text(size = 12), strip.text.x= element_text(size=12))+xlab("Efficacy level") + ylab("Amplitude [µV]") + #ggtitle("Late CNV")+  #ylim(-0.12, 0.15)+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.position="none")

pdf("Figures/ERN_Efficacy_Reward_Acc.pdf", width =4, height = 4)
pRT
dev.off()

eff_df <- Effect(c("cRT", "Congruency"), ERNmod0Con, xlevels=list(cRT =seq(min(a2$cRT, na.rm=TRUE ), max(a2$cRT, na.rm=TRUE), 0.01)))
IA <- as.data.frame(eff_df)
IA$Congruency <- ordered(IA$Congruency, levels=c( "Incongruent", "Neutral", "Congruent"))
pRT<- ggplot(data=IA, aes(x=cRT, y=fit, color=Congruency  )) + geom_line()+scale_colour_manual(name="Congruency", values=cbPalette) +theme_bw(12)+ geom_ribbon(data=IA, aes(x=cRT, max = fit + se, min = fit- se, fill = Congruency),alpha=0.1, inherit.aes = FALSE)+
  scale_fill_manual(name="Congruency", values=cbPalette) + xlab("RT") + ylab("ERN") + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) #+ theme(legend.position=c(0.25, 0.2))

pdf("Figures/ERN_Congruency_RT.pdf", width =4, height = 4)
pRT
dev.off()

## errors only
print (summary(ERNmod0Conerr <- lmer(ERN~ (EffLvl*RewLvl)+Congruency*((cRT +I(cRT*cRT)))+(meanCorr1) +crBaseline+(  I(cRT*cRT)|SubID), a2[!a2$RT<200 & a2$Acc==0,], 
                                  REML=FALSE)))

sv1_max <- svd(getME(ERNmod0Conerr, "Tlist")[[1]]) 
sv1_max$d
round(sv1_max$d^2/sum(sv1_max$d^2)*100, 1)  

print (summary(ERNmod0Concorr <- lmer(ERN~ (EffLvl*RewLvl)+Congruency*((cRT +I(cRT*cRT)))+(meanCorr1) +crBaseline+(  I(cRT*cRT)|SubID), a2[!a2$RT<200 & a2$Acc==1,], 
                                     REML=FALSE)))

sv1_max <- svd(getME(ERNmod0Concorr, "Tlist")[[1]]) 
sv1_max$d
round(sv1_max$d^2/sum(sv1_max$d^2)*100, 1)  

tab_model(ERNmod0Conerr, ERNmod0Concorr)


# 5b) FRN

print (summary(FRNmod0CP2P <- lmer(P2PFRN~ IsRewarded*(EffLvl+RewLvl)+ Congruency+cTrial+(IsRewarded|SubID), a1[!a1$IsMiss & !a1$RT<200,], 
                                REML=FALSE)))

sv1_max <- svd(getME(FRNmod0CP2P, "Tlist")[[1]]) 
sv1_max$d
round(sv1_max$d^2/sum(sv1_max$d^2)*100, 1) 

tab_model(FRNmod0CP2P)

sum1 <-summarySEwithin(a1[!a1$IsMiss & !a1$RT<200,], measurevar="P2PFRN", withinvars=c("EffLvl", "RewLvl", "IsRewarded"),  idvar="SubID", na.rm=TRUE)

pFRNEff <- ggplot(data=sum1, aes(x=EffLvl, y=P2PFRN , group = interaction(IsRewarded, RewLvl), color = IsRewarded,  fill= RewLvl)) +  geom_errorbar(aes(max = P2PFRN+se, min = P2PFRN-se), width=0.1,position=position_dodge(.1)) +geom_point(shape = 22, size = 5,position=position_dodge(.1))+geom_line(position=position_dodge(.1))+ scale_fill_manual(name = "Reward level", values=c("#343085" ,"#F1E923"))+
  theme_bw()+ scale_colour_manual(name="Is Rewarded",values=c("#CC0000","#000000"))+#facet_wrap(~EffLvl)+
  theme(axis.text = element_text(size = 12), strip.text.x= element_text(size=12))+xlab("Efficacy level") + ylab("Amplitude [µV]") + #ggtitle("Late CNV")+  #ylim(-0.12, 0.15)+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) #+ theme(legend.position="none")


pdf("Figures/FRN_Reward_and_Efficacy_by_reward_within_SE_all_in_one.pdf", width =5, height = 4)
pFRNEff
dev.off()

## get residuals from model without trial and plot them as function of trial (Supplement)
print (summary(AccRTmod0tbnt <- lmer(AccRT~ EffLvl*RewLvl +Congruency+(Congruency|SubID), a1[!a1$IsMiss & !a1$RT<200,], 
                                     REML=FALSE)))

sv1_max <- svd(getME(AccRTmod0tbnt, "Tlist")[[1]]) 
sv1_max$d
round(sv1_max$d^2/sum(sv1_max$d^2)*100, 1)  

a1$AccRTResNT <- rep(NA_real_, length(a1$SubID))

a1$AccRTResNT[!a1$IsMiss & !a1$RT<200 & !is.nan(a1$AccRT)] <- resid(AccRTmod0tbnt)

sum1 <-summarySEwithin(a1[!a1$IsMiss & !a1$RT<200,], measurevar="AccRTResNT", withinvars=c("Trial"),  idvar="SubID", na.rm=TRUE)

sum1$Trial <- as.numeric(sum1$Trial)
pDRmCNVl <- ggplot(data=sum1, aes(x=Trial, y=AccRTResNT))+geom_ribbon(aes(max = AccRTResNT+se, min = AccRTResNT-se)) +  geom_point(shape = 22, size = 1)+ scale_fill_manual(name = "Reward level", values=c("#343085" ,"#F1E923"))+
  theme_bw()+ geom_hline(aes(yintercept=0), colour="#000000", linetype=2, size=0.4)+ geom_smooth()+ #facet_wrap(~ Congruency)+
  theme(axis.text = element_text(size = 12), strip.text.x= element_text(size=12))+xlab("Trial") + ylab("Average Residual RT") + #ggtitle("Late CNV")+  #ylim(-0.12, 0.15)+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.position="none")

pdf("Figures/ResidualRTTrial.pdf", width =4, height = 4)
pDRmCNVl
dev.off()

### look at topos of effects

ap <- a1

ap  <- cbind(a1,P3b )

plotEmat5 <- as.data.frame(matrix(, nrow = 35, ncol=4))
plotTmat5 <- as.data.frame(matrix(, nrow = 35, ncol=4))

j=1
for (i in (ncol(a1)+1):(ncol(ap))){
  ap$DV  <- ap[,i]
  Plotmod <- lmer(DV~ EffLvl*RewLvl+(1|SubID),ap[!ap$IsMiss==1 & !ap$RT<200,],
                  REML=FALSE)
  k <- summary(Plotmod)
  if (i ==(ncol(a1)+1)){
    f <- dimnames(k$coefficients)[[1]]
    colnames(plotEmat5) <-  f
    colnames(plotTmat5) <-  f}
  plotEmat5[j,] <- k$coefficients[,1]
  plotTmat5[j,] <- k$coefficients[,3]
  j=j+1
}

colnames(plotEmat5) <- sedit(colnames(plotEmat5), "EffLvl2-1", "Efficacy")
colnames(plotEmat5) <- sedit(colnames(plotEmat5), "RewLvl2-1", "Reward")
colnames(plotEmat5) <- sedit(colnames(plotEmat5), ":", "by")
colnames(plotTmat5) <- sedit(colnames(plotTmat5), "EffLvl2-1", "Efficacy")
colnames(plotTmat5) <- sedit(colnames(plotTmat5), "RewLvl2-1", "Reward")
colnames(plotTmat5) <- sedit(colnames(plotTmat5), ":", "by")


writeMat('Export/plotes5.mat',Emat5=plotEmat5)
writeMat('Export/plotts5.mat',Tmat5=plotTmat5)

## CNV
ap <- a1

ap  <- cbind(a1, CNV)
plotEmat2 <- as.data.frame(matrix(, nrow = 35, ncol=4))
plotTmat2 <- as.data.frame(matrix(, nrow = 35, ncol=4))

j=1
for (i in (ncol(a1)+1):(ncol(ap))){
  ap$DV  <- ap[,i]
  Plotmod <- lmer(DV~ EffLvl*RewLvl+(1|SubID),ap[!ap$IsMiss==1 & !ap$RT<200,], 
                  REML=FALSE)
  k <- summary(Plotmod)
  if (i ==(ncol(a1)+1)){
    f <- dimnames(k$coefficients)[[1]]
    colnames(plotEmat2) <-  f
    colnames(plotTmat2) <-  f}
  plotEmat2[j,] <- k$coefficients[,1]
  plotTmat2[j,] <- k$coefficients[,3]
  j=j+1
}

colnames(plotEmat2) <- sedit(colnames(plotEmat2), "EffLvl2-1", "Efficacy")
colnames(plotEmat2) <- sedit(colnames(plotEmat2), "RewLvl2-1", "Reward")
colnames(plotEmat2) <- sedit(colnames(plotEmat2), ":", "by")
colnames(plotTmat2) <- sedit(colnames(plotTmat2), "EffLvl2-1", "Efficacy")
colnames(plotTmat2) <- sedit(colnames(plotTmat2), "RewLvl2-1", "Reward")
colnames(plotTmat2) <- sedit(colnames(plotTmat2), ":", "by")


writeMat('Export/plotes2.mat',Emat2=plotEmat2)
writeMat('Export/plotts2.mat',Tmat2=plotTmat2)

######################## Theta r-locked #####################################

# load and merge subject files
logdir= 'Timefreq/'

fl <- list.files(logdir)

for (i in 1:length (fl))
{ 

tmp <- read.csv(paste0(logdir, fl[i]))#, header = TRUE, sep = "", #dec = ",", bad... kills everything! fill = TRUE)
if (i==1)
{a3 <- tmp}
else {a3 <- rbind(a3,tmp)}
}

## exclude bad subs
a3 <- a3[!a3$subid==35 & !a3$subid==1 & !a3$subid==21 & !a3$subid==38 & !a3$subid==39 & !a3$subid==43 & !a3$subid==46 & !a3$subid==53 & !a3$subid==63,]

## prep variables

a3$subid <- as.factor(a3$subid)
a3$ctrial <- scale(a3$trial, scale= FALSE, center = TRUE)
a3$fAcc <- as.factor(a3$acc)
contrasts(a3$fAcc) <- contr.sdif(2)
a3$RewLvl <- as.factor(a3$rewlvl)
a3$RewLvl <- revalue(a3$RewLvl , c("-0.5"="1", "0.5"="2"))
contrasts(a3$RewLvl) <- contr.sdif(2)

a3$EffLvl <- as.factor(a3$efflvl)
a3$EffLvl <- revalue(a3$EffLvl , c("-0.5"="0", "0.5"="1"))
contrasts(a3$EffLvl) <- contr.sdif(2)

a3$Congruency <- as.factor(a3$congruency)
a3$Congruency <- revalue(a3$Congruency , c("-0.5"="Incongruent", "0.5"="Congruent", "0"="Neutral"))
a3$Congruency<- ordered(a3$Congruency, levels=c("Incongruent", "Neutral","Congruent"))
contrasts(a3$Congruency) <- contr.sdif(3)

a3$cfcz_resp_theta_baseline <- scale(a3$fcz_resp_theta_baseline, scale=TRUE, center=TRUE)

a3$crt <- scale(a3$rt, scale = FALSE, center = TRUE)

print (summary(Thetamod0Con <- lmer(fcz_resp_theta~ (EffLvl+RewLvl)*fAcc+cfcz_resp_theta_baseline+(fAcc |subid), a3[!a3$rt<200,], 
                                  REML=FALSE)))

sv1_max <- svd(getME(Thetamod0Con, "Tlist")[[1]]) 
sv1_max$d
round(sv1_max$d^2/sum(sv1_max$d^2)*100, 1)  


tab_model(Thetamod0Con)


print (summary(Thetamod0Connest <- lmer(fcz_resp_theta~ fAcc/(EffLvl+RewLvl)+cfcz_resp_theta_baseline+(fAcc |subid), a3[!a3$rt<200,], 
                                    REML=FALSE)))

sv1_max <- svd(getME(Thetamod0Connest, "Tlist")[[1]]) 
sv1_max$d
round(sv1_max$d^2/sum(sv1_max$d^2)*100, 1)  


tab_model(Thetamod0Connest)



##### Between group comparison #####

mergedData <- rbind(a1[,c(1:8,10:22)], a1a[,1:21])

mergedData$Vers <- as.factor(mergedData$Vers)
contrasts(mergedData$Vers) <- contr.sdif(2)

mergedData$SubID <- as.factor(mergedData$SubID)

mergedData$cTrial <- scale(mergedData$Trial/100, scale= FALSE, center=TRUE)

mergedData$oRewLvL <- mergedData$RewLvl

mergedData$RewLvl[mergedData$RewLvl==4] <- 2

mergedData$RewLvl <- factor(mergedData$RewLvl)
contrasts(mergedData$RewLvl) <- contr.sdif(2)
contrasts(mergedData$EffLvl) <- contr.sdif(2)
mergedData$Congruency<- ordered(mergedData$Congruency, levels=c("Incongruent", "Neutral","Congruent"))

contrasts(mergedData$Congruency) <- contr.sdif(3)

print (summary(ACCmod0 <- glmer(Acc~ (EffLvl*RewLvl + Congruency + cTrial)*Vers+(Congruency|SubID), data = mergedData[!mergedData$IsMiss & !mergedData$RT<200,], 
                                family = binomial))) 

sv1_max <- svd(getME(ACCmod0, "Tlist")[[1]]) 
sv1_max$d
round(sv1_max$d^2/sum(sv1_max$d^2)*100, 1)  


print (summary(AccRTmod0t <- lmer(AccRT~ (EffLvl*RewLvl +Congruency+cTrial)*Vers+(EffLvl+Congruency|SubID), mergedData[!mergedData$IsMiss & !mergedData$RT<200,], 
                                  REML=FALSE)))

sv1_max <- svd(getME(AccRTmod0t, "Tlist")[[1]]) 
sv1_max$d
round(sv1_max$d^2/sum(sv1_max$d^2)*100, 1)  

tab_model(ACCmod0, AccRTmod0t,transform=NULL )


