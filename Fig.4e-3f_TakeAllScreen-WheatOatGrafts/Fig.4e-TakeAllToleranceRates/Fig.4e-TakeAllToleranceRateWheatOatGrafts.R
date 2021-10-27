##############################################
##############################################
# R-script by Gregory Reeves (27 MAY 2021)
# University of Cambridge, Department of Plant Sciences,
# Downing Site, Cambridge CB2 3EA, United Kingdom
# E-mail me at drgregreeves@gmail.com or gr360@cam.ac.uk
##############################################
##############################################



library(tidyverse)
library(ggpubr)
library(dplyr)



#######################
# Plot boxplot for the resistance data
######################


### You need to specify the working directory. 
### You can always use the command "file.choose()", then copy the working directory with the file "DiseaseScreenResults.csv"
setwd("~/")

DiseaseScreenResults <- read.csv("DiseaseScreenResults.csv")

### Alternatively use:
# file.choose()


#######################
AnBarInDen<- ggplot(DiseaseScreenResults, 
                    aes(x=GraftType, 
                        y=ToleranceRate, 
                        fill=GraftType,
                        linetype="solid"))+
  geom_bar(position=position_dodge(), 
           stat="identity",
           width=0.75) +
  scale_y_continuous(name=expression("Take All Tolerance Rate (%)"), expand=c(0,0)) + 
  scale_fill_manual(values=c("#2ca089ff", "#677821ff", "#502d16ff")) +
  theme(axis.line = element_line(linetype = "solid"), 
        axis.ticks = element_line(), 
        axis.ticks.length = unit(-0.25 , "cm"), # -ve length = inside ticks
        panel.grid.major = element_line(colour = NA), 
        panel.grid.minor = element_line(colour = NA), 
        axis.title = element_text(size = 18,
                                  face = "plain"), 
        axis.text = element_text(size = 15, 
                                 face = "plain", 
                                 colour = "black"), 
        axis.text.x = element_text(size = 15,
                                   margin=margin(10,5,10,5,"pt")), 
        axis.text.y = element_text(size = 15,
                                   margin=margin(5,10,10,5,"pt")), 
        legend.text = element_text(size = 9), 
        legend.title = element_text(size = 10, 
                                    face = "bold", 
                                    colour = "black"), 
        panel.background = element_rect(fill = NA), 
        legend.position="none", #Un-comment to remove the legend
        legend.key = element_rect(fill = NA), 
        legend.background = element_rect(fill = NA), 
        legend.direction = "horizontal") +
  labs(x = "Graft Type")
AnBarInDen


