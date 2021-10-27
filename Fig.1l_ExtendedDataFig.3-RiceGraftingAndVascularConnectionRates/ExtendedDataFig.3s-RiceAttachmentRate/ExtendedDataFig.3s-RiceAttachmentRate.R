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


#############################################
# Attachment Rates 
#############################################

### You need to specify the working directory. 
### You can always use the command "file.choose()", then copy the working directory with the file "AttachmentData.csv"
#setwd("~/")


file.choose() ### Choose file called "AttachmentData.csv"


#Grafting data input
AttachmentData <- read.csv("AttachmentData.csv", sep=",",header=TRUE)


# Rename Columne in AttachmentData
colnames(AttachmentData) <-c("Days_after_grafting",	
                             "Type",	
                             "Sample",	
                             "Not_attached",
                             "Attached",
                             "number",	
                             "AttachmentRate")


# Summarise data
AttachmentDataSummary <- summarise(group_by(AttachmentData, 
                                            Days_after_grafting, Type),
                                   meanAttachmentRate=mean(AttachmentRate),
                                   sdAttachmentRate=sd(AttachmentRate),
                                   seAttachmentRate=sd(AttachmentRate)/sqrt(n()),
                                   SumNumber=sum(number),
                                   n=n())


# Code to plot the data
AtachmentPlot <- ggplot(data=AttachmentDataSummary,
                        mapping=aes(Days_after_grafting,
                                    meanAttachmentRate,
                                    group=Type, 
                                    color=Type)) +
  #geom_line() + 
  geom_errorbar(aes(Days_after_grafting, 
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
  scale_x_continuous(breaks=seq(0,14,1))+
  scale_y_continuous(breaks=seq(0, 100, 5))+
  labs(x = expression(Days~after~grafting),
       y = expression(Attachment~rate~("%"))) +
  ### Generates the line of best fit and inlcudes the Standard Error shading  THIS IS NOT THE BINOMIAL REGRESSION LINE    
  #geom_smooth(se = FALSE, size = 1.5 #method = "auto" formula = y ~ log(x))+
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
AtachmentPlot






