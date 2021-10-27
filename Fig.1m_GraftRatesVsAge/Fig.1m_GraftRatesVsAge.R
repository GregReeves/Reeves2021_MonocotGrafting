###################################################
# Grafting Rates for Rice Self Grafts Over Time
# R-script by Greg Reeves (16 June 2020, Would be Dad's birthday)
###################################################

library(tidyverse)
library(ggplot2)
library(plotly)
library(mgcv) 
library(dplyr)
library(ggpubr)
library(ggrepel)
library(ggThemeAssist)
library(car)
library(multcompView)
library(rstatix)
library(lsmeans)

#############################################
# Attachment Rates Over Time
#############################################



### You need to specify the working directory. 
### You can always use the command "file.choose()", then copy the working directory with the file "RiceTimecourseTranscriptome-21Oct2020.csv"
#setwd("~/")

filelocation <- file.choose() ### Choose the file called "TimeCourse-RiceGrafts14June2020.csv"


#Grafting data input
TimeCourseData <- read.csv(filelocation, sep=",",header=TRUE)




#####################
#We need to filter out the controls, these mess up the calculation

# Rename Columne in AttachmentData
colnames(TimeCourseData) <-c("DaysAfterGermination",	
                             "Type",
                             "Sample",
                             "Attached",
                             "NoAttached",	
                             "Total",
                             "AttachmentRate")

#####################
#We need to filter out the controls, these mess up the calculation
TimeCourseFILTERED <- filter(TimeCourseData, !(Type == 'Control'))

#####################
#  library("ggpubr")
ggline(TimeCourseFILTERED , x = "DaysAfterGermination", y = "AttachmentRate", 
       add = c("mean_sd", "jitter"), 
       #order = c("ctrl", "trt1", "trt2"),
       ylab = "Attachment Rate (%)", xlab = "Days After Germination")
#####################
#####################



# Summarise data
TimeCourseFILTEREDSummary <- summarise(group_by(TimeCourseFILTERED, 
                                                DaysAfterGermination, Type),
                                       meanAttachmentRate=mean(AttachmentRate),
                                       sdAttachmentRate=sd(AttachmentRate),
                                       seAttachmentRate=sd(AttachmentRate)/sqrt(n()),
                                       n=n())

####################
AtachmentTimecoureFILTEREDPlot <- ggplot(data=TimeCourseFILTEREDSummary,
                                         mapping=aes(DaysAfterGermination,
                                                     meanAttachmentRate,
                                                     group=Type, 
                                                     color=Type)) +
  geom_errorbar(aes(DaysAfterGermination, 
                    ymin=meanAttachmentRate-sdAttachmentRate, 
                    ymax=meanAttachmentRate+sdAttachmentRate,
                    width=0.2)) +
  geom_point(aes(shape=Type, 
                 color=Type, 
                 fill=Type), 
             size=4) + 
  scale_shape_manual(values=c(21))+ 
  scale_fill_manual(breaks = c("Graft"),
                    values=c("#f33829ff"))+
  #Places a tick every 2 days from 0 to 10 days
  scale_x_continuous(breaks=seq(0,10,1))+
  scale_y_continuous(breaks=seq(0, 60, 10))+
  labs(x = expression(Days~after~germination),
       y = expression(Attachment~rate~("%"))) +
  # Generates the line of best fit and inlcudes the Standard Error shading
  geom_smooth(se = TRUE,
              size = 1.5,
              method = "auto"
              #formula = y ~ log(x)
  )+
  scale_color_manual(breaks = c("Graft"), values=c("#ff7751ff"))+
  theme_classic()+
  theme(legend.position=c(0.75, 0.8),
        axis.line.x = element_line(colour = "#000000", size = 1.5, linetype = "solid"),
        axis.line.y = element_line(colour = "#000000", size = 1.5, linetype = "solid"),
        axis.ticks = element_line(colour = "#000000", size = 1.5, linetype = "solid"),
        axis.text.x = element_text(face="bold", size=25, colour = "#333333"),
        axis.text.y = element_text(face="bold", size=25, colour = "#333333"),
        axis.title = element_text(size=30)) 

# Display AtachmentPlot with Classic settings
AtachmentTimecoureFILTEREDPlot


###############
# Binomial regression on the data
##############


library(lme4)
#library(ggThemeAssist) #ggThemeAssistGadget(plot)

#Get Data
dat<-read.csv(file.choose()) ### Choose the file called "TimeCourse-RiceGrafts14June2020.csv"



#####################
#We need to filter out the controls, these mess up the calculation
datFILTERED <- filter(dat, !(Type == 'Control'))

#Get Data
dat<-datFILTERED

#Add rename same variables

dat$Days<-dat$DaysAfterGerm

#Fit binomial mixed effects model with Days as a fixed effect predictor and Sample as a random effect
mod<-glmer(cbind(Attached , NotAttached)~Days+(1|Sample) , data=dat , family='binomial')

#Fit binomial mixed effects model only Sample as a random effect
nullmod<-glmer(cbind(Attached , NotAttached)~(1|Sample) , data=dat , family='binomial')

#Compare the two models with a Likelihood ratio test
anova(mod , nullmod , test="LRT") #### THIS IS THE P-VALUE IN Fig. 1m

#The estimated constant attachment rate can be calculated from the intercept coefficient in the following table
summary(nullmod)

#This helps extract the log(odds) of the attachment probability
summaryobj <- summary(nullmod)
coeff <- summaryobj$coefficients
coeff[1] #This is the Estimate value
coeff[2] #This is the standard error
est<-family(nullmod)$linkinv(coeff[1])
est
#So the constant attachment rate is about 0.45. We can also calucate the upper and lower 95% confidence intervals using the Std.Error value (0.07545)
#from the summary table:

upper<-family(nullmod)$linkinv(coeff[1] + 1.96*coeff[2])
lower<-family(nullmod)$linkinv(coeff[1] - 1.96*coeff[2])


plot(jitter(dat$Days) , dat$Attached/dat$Total)
abline(h=est , col='black')
abline(h=c(lower,upper) , col='red' , lty=2)


#create the basic glm without random effects
mod<-glm(cbind(Attached , NotAttached)~Days , data=dat , family='binomial')

#create a set of points to estimate the curves at
predDat<-data.frame(Days=1:10)

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



#######################
# Let's see a pretty plot with the non-grafted control
#######################



# Summarise data
TimeCourseDataSummary <- summarise(group_by(TimeCourseData, 
                                            DaysAfterGermination, Type),
                                   meanAttachmentRate=mean(AttachmentRate),
                                   sdAttachmentRate=sd(AttachmentRate),
                                   seAttachmentRate=sd(AttachmentRate)/sqrt(n()),
                                   n=n())


#######################
# Code to plot the data
#######################

AtachmentTimecourePlot <- ggplot(data=TimeCourseDataSummary,
                                 mapping=aes(DaysAfterGermination,
                                             meanAttachmentRate,
                                             group=Type, 
                                             color=Type)) +
  #geom_line() + 
  geom_errorbar(aes(DaysAfterGermination, 
                    ymin=meanAttachmentRate-sdAttachmentRate, 
                    ymax=meanAttachmentRate+sdAttachmentRate,
                    width=0.2)) +
  geom_point(aes(shape=Type, 
                 color=Type, 
                 fill=Type), 
             size=4) + 
  scale_shape_manual(values=c(22,21))+ 
  scale_fill_manual(breaks = c("Control", "Graft"),
                    values=c("#0eb9cbff", "#f33829ff"))+
  #Places a tick every 2 days from 0 to 10 days
  scale_x_continuous(breaks=seq(0,10,1))+
  scale_y_continuous(breaks=seq(0, 100, 10))+
  labs(x = expression(Days~after~germination),
       y = expression(Attachment~rate~("%"))) +
  # Generates the line of best fit and inlcudes the Standard Error shading
  geom_smooth(se = TRUE,
              size = 1.5,
              method = "auto"
              #formula = y ~ log(x)
  )+
  scale_color_manual(breaks = c("Control", "Graft"), values=c("#027381ff", "#ff7751ff"))+
  theme_classic()+
  theme(legend.position=c(0.75, 0.8),
        axis.line.x = element_line(colour = "#000000", size = 1.5, linetype = "solid"),
        axis.line.y = element_line(colour = "#000000", size = 1.5, linetype = "solid"),
        axis.ticks = element_line(colour = "#000000", size = 1.5, linetype = "solid"),
        axis.text.x = element_text(face="bold", size=25, colour = "#333333"),
        axis.text.y = element_text(face="bold", size=25, colour = "#333333"),
        axis.title = element_text(size=30)) 

# Display AtachmentPlot with Classic settings
AtachmentTimecourePlot
