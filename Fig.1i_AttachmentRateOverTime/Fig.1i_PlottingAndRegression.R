##############################################
##############################################
# R-script by Gregory Reeves (27 MAY 2021)
# University of Cambridge, Department of Plant Sciences,
# Downing Site, Cambridge CB2 3EA, United Kingdom
# E-mail me at drgregreeves@gmail.com or gr360@cam.ac.uk
##############################################
##############################################
library(tidyverse)
library(ggplot2)
library(mgcv) 
library(dplyr)
library(ggpubr)
library(car)
library(rstatix)
library(lsmeans)
library(lme4) # You need this package for the binomial regression steps

###################################################
# Grafting Rates for Rice Self Grafts Over Time
###################################################

### You need to specify the working directory. 
### You can always use the command "file.choose()", then copy the working directory with the file "RiceTimecourseTranscriptome-21Oct2020.csv"
setwd("~/")

file.choose()


#Grafting data input
TranscriptomeData <- read.csv("RiceTimecourseTranscriptome-21Oct2020.csv", sep=",",header=TRUE)

colnames(TranscriptomeData) <-c("DateGrafted",
                                "DateEvaluated",
                                "DaysAfterGrafting",
                                "DaysAfterGermination",
                                "Sample",
                                "RiceGenotype",
                                "GraftType",
                                "Fused",
                                "NotSeparated",
                                "Separated",
                                "Total",
                                "AttachmentEfficiency",
                                "SpecialNotes")


#####################
#We need to filter out the controls, these mess up the calculation
TranscriptomeDataFILTERED <- filter(TranscriptomeData, (GraftType == 'SelfGraft'))
TranscriptomeDataFILTERED2 <- filter(TranscriptomeDataFILTERED, !(DaysAfterGermination == '10'))

#####################
#  library("ggpubr")
AttachmentPlot <-ggline(TranscriptomeDataFILTERED2, 
                        x = "DaysAfterGrafting", 
                        y = "AttachmentEfficiency", 
                        add = c("mean_se"), 
                        add.params = list(color = "#e56814", size= 0.5),
                        color = c("#c45911"), 
                        point.color = c("#a34a0e"),
                        #order = c("ctrl", "trt1", "trt2"),
                        ylab = "Attachment Rate (%)", xlab = "Days After Grafting",
                        size = 1) +
  scale_x_continuous(breaks=seq(1,7,by=1)) +
  scale_y_continuous(breaks=seq(0,80,by=10))
#####################
AttachmentPlot + theme(axis.line = element_line(size = 2.5, 
                                                linetype = "solid"), 
                       axis.ticks = element_line(size = 2.5), 
                       panel.grid.major = element_line(colour = NA), 
                       panel.grid.minor = element_line(colour = NA), 
                       axis.title = element_text(size = 30, 
                                                 face = "bold"), 
                       axis.text = element_text(size = 30, 
                                                face = "bold", colour = "gray30"), 
                       axis.text.x = element_text(size = 27), 
                       axis.text.y = element_text(size = 27), 
                       legend.text = element_text(size = 27), 
                       legend.title = element_text(size = 27, 
                                                   face = "bold", colour = "gray30"), 
                       panel.background = element_rect(fill = NA), 
                       legend.key = element_rect(fill = NA), 
                       legend.background = element_rect(fill = NA))
#####################


###############
# Mixed Effect Binomial regression on the data
##############


# library(lme4)

### You need to specify the working directory. 
### You can always use the command "file.choose()", then copy the working directory with the file "RiceTimecourseTranscriptome-21Oct2020.csv"
setwd("~/")




#Grafting data input
dat <- read.csv("RiceTimecourseTranscriptome-21Oct2020.csv", sep=",",header=TRUE)

colnames(dat) <-c("DateGrafted",
                  "DateEvaluated",
                  "DaysAfterGrafting",
                  "DaysAfterGermination",
                  "Sample",
                  "RiceGenotype",
                  "Type",
                  "Attached",
                  "NotSeparated",
                  "Separated",
                  "Total",
                  "AttachmentEfficiency",
                  "SpecialNotes")


#Add rename same variables

dat$Days<-dat$DaysAfterGrafting
dat$NotAttached<-dat$Total-dat$Attached

#####################
#We need to filter out the controls, these mess up the calculation
datFILTERED <- filter(dat, !(Type == 'NonGraftControl'))
datFILTERED2 <- filter(datFILTERED, !(DaysAfterGermination == '10'))

#Get Data
dat<-datFILTERED2



#Fit binomial mixed effects model with Days as a fixed effect predictor and Sample as a random effect
mod<-glmer(cbind(Attached , NotAttached)~Days+(1|Sample) , data=dat , family='binomial')

#Fit binomial mixed effects model only Sample as a random effect
nullmod<-glmer(cbind(Attached , NotAttached)~(1|Sample) , data=dat , family='binomial')

#Compare the two models with a Likelihood ratio test
anova(mod , nullmod , test="LRT") # This yields the P-value for the plot in Fig. 1i


#The estimated constant attachment rate can be calculated from the intercept coefficient in the following table
summary(nullmod)

#This helps extract the log(odds) of the attachment probability
summaryobj <- summary(nullmod)
coeff <- summaryobj$coefficients
coeff[1] #This is the Estimate value
coeff[2] #This is the standard error
#This value (Estimate) isn't the attachment probability directly but is instead the log(odds) of the attachment
#probability. We can recover the exact value through the following function

est<-family(nullmod)$linkinv(coeff[1])
est
#We can also calculate the upper and lower 95% confidence intervals using the Std.Error value
#from the summary table:

upper<-family(nullmod)$linkinv(coeff[1] + 1.96*coeff[2])
lower<-family(nullmod)$linkinv(coeff[1] - 1.96*coeff[2])


#The basic visualisation for this null model ("null" because there are no statistically significant fixed effects) is below/ I've added a small
#amount of jitter to the x-component of each point to prevent overplotting on some of the days (you have multiple measurements that are exactly
#the same on some days and this afct would be obscured if you just plotted the data directly).

plot(jitter(dat$Days) , dat$Attached/dat$Total)
abline(h=est , col='black')
abline(h=c(lower,upper) , col='red' , lty=2)

#FYI the visualisation for the model with Days as a predictor variable is below, just so you have the
#code. As I mentioned before, I can only add the 95% confidence lines here because there isn't any random effect due to Sample. Ordinarily you'd just
#be stuck with the fitted line.

#create the basic glm without random effects
mod<-glm(cbind(Attached , NotAttached)~Days , data=dat , family='binomial')

#create a set of points to estimate the curves at
predDat<-data.frame(Days=1:7)

#Predictions and std errors on the linearised log-odds scale
est<-predict(mod , predDat , se.fit=TRUE)$fit
stderr<-predict(mod , predDat , se.fit=TRUE)$se.fit

#Predictions and 95% confidence intervals on the probability scale
meanPred<-family(mod)$linkinv(est)
upperPred<-family(mod)$linkinv(est + 1.96*stderr)
lowerPred<-family(mod)$linkinv(est - 1.96*stderr)


plot(jitter(dat$Days) , dat$Attached/dat$Total)+
  points(meanPred~predDat$Days , type='l' , col='black')
points(upperPred~predDat$Days , type='l' , col='red' , lty=2)
points(lowerPred~predDat$Days , type='l' , col='red' , lty=2)
