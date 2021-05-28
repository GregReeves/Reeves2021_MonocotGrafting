##############################################
##############################################
# R-script by Gregory Reeves (27 MAY 2021)
# University of Cambridge, Department of Plant Sciences,
# Downing Site, Cambridge CB2 3EA, United Kingdom
# E-mail me at drgregreeves@gmail.com or gr360@cam.ac.uk
##############################################
##############################################


library(tidyverse)
library(viridis)
library(ggrepel)
library(ggpubr)
# library(ggplot2)



### You need to specify the working directory. 
### You can always use the command "file.choose()", then copy the working directory with the file "HormoneRiceGraftBoost_Dec2020.csv"
setwd("~/")



#Grafting data input
HORMONEINHIBORData <- read.csv("HormoneRiceGraftBoost_Dec2020.csv", sep=",",header=TRUE)

HORMONEINHIBORData$Treatment <- as.character(HORMONEINHIBORData$Treatment)

HORMONEINHIBORDataSummary <- HORMONEINHIBORData  %>%
  group_by(Treatment) %>%
  summarise(meanGraftingRate=100*(sum(No.Fused)/sum(No.Germinated)), 
            sumNoFused=sum(No.Fused),
            sumNoTotal=sum(No.Germinated), 
            FusionRate = FusionRate,
            TreatmentType =TreatmentType,
            sdGraftingRate=sd(FusionRate))

#####################


ManualOrder <- c("HalfMS_Hormone", 
                 "HalfMS_1mgLGA3", "HalfMS_1uM2_4D","HalfMS_1uM2_4D_1mgLGA3",
                 "HalfMS_DMSO",
                 "HalfMS_100uM_PBZ", "HalfMS_100uM_TIBA", "HalfMS_100uM_TIBA_PBZ")



ggboxplot(data = HORMONEINHIBORData, 
          x="Treatment", 
          y="FusionRate", 
          #order=AnatOrderList, 
          fill= "Treatment", 
          color = "black", 
          #add = "jitter",
          weight="5", 
          order = ManualOrder,
          #add = "jitter",  
          add.params = list(color = "Treatment")) +
  scale_y_continuous(name= "Grfating rate (%)", 
                     breaks = seq(0,100,by=5),
                     limits = c(0, 100)) +
  scale_fill_viridis_d() +
  facet_grid(~ TreatmentType, scales = "free_x", space = "free_x") + # Comment this out to intercalate the bars
  ###Minimal Theme
  theme(axis.line = element_line(linetype = "solid"), 
        axis.ticks = element_line(), 
        #minor.tick(nx=2, ny=2, tick.ratio=0.5, x.args = list(), y.args = list()),
        axis.ticks.length = unit(-0.25 , "cm"), # -ve length = inside ticks
        panel.grid.major = element_line(colour = NA), 
        panel.grid.minor = element_line(colour = NA), 
        axis.title = element_text(size = 18,
                                  face = "plain"), 
        axis.text = element_text(size = 15, 
                                 face = "plain", 
                                 colour = "black"), 
        axis.text.x = element_text(size = 15,
                                   margin=margin(10,5,10,5,"pt"),
                                   angle=45), 
        axis.text.y = element_text(size = 15,
                                   margin=margin(5,10,10,5,"pt")), 
        legend.text = element_text(size = 9), 
        legend.title = element_text(size = 15, 
                                    face = "bold", 
                                    colour = "black"), 
        panel.background = element_rect(fill = NA), 
        legend.position="none",
        legend.key = element_rect(fill = NA), 
        legend.background = element_rect(fill = NA), 
        legend.direction = "horizontal") +labs(x = "Treatment", 
                                               y = "Grafting rate (%)")



