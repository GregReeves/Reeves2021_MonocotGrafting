##############################################
##############################################
# R-script by Gregory Reeves (27 MAY 2021)
# University of Cambridge, Department of Plant Sciences,
# Downing Site, Cambridge CB2 3EA, United Kingdom
# E-mail me at drgregreeves@gmail.com or gr360@cam.ac.uk
##############################################
##############################################


#######################
# Plots for grafting rates among various monocot species
# And plots for grafts between C3 and C4 species
#######################

library(ggplot2)
library(tidyverse)
library(viridis)
#library(Hmisc) #Allows minor.ticks function in ggplot



### You need to specify the working directory. 
### You can always use the command "file.choose()", then copy the working directory with the file "Interspecific C3-C4 grafts.csv"
setwd("~/")



#Grafting data input
C3C4graftdata <- read.csv("Interspecific C3-C4 grafts.csv", sep=",",header=TRUE)

#Combines the SCION and ROOTSTOCK into new column
C3C4graftdata$SpecificCombo <- paste0(C3C4graftdata$SCION, 
                                      C3C4graftdata$ROOTSTOCK)

C3C4graftdata_summarised <- C3C4graftdata %>%
  group_by(Type, 
           SpecificCombo) %>%
  summarise(meanPercentFused=mean(PercentFused),
            Total_n=sum(No.GraftAttempts))


#####################
#We need to filter by combo to make each plot
### ADD WHICHEVER SPECIES COMBINATION YOU WANT HERE
SpecificComboFILTERED <-  filter(C3C4graftdata_summarised, (Type == 'Wheat-DurumWheat'))



#####################
# Make the plot 

ggplot(data = SpecificComboFILTERED, aes(x=SpecificCombo, 
                                         y=meanPercentFused, 
                                         fill = SpecificCombo)) + 
  geom_bar(stat="identity", position="identity", width = 0.9) +
  geom_text(aes(label = paste("n =", Total_n, sep=" "), 
                vjust = -1),
            color = "gray30",
            size = 4) + 
  scale_y_continuous(breaks=seq(0,50,by=5),
                     limits=c(0,50)) +
  #scale_x_continuous(breaks=seq(1,7,by=2)) +
  scale_fill_manual(values = c("#999999", "#E69F00")) +
  #scale_fill_viridis_d() + 
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
                                   margin=margin(10,5,10,5,"pt")), 
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
        legend.direction = "horizontal") +labs(x = "Species combination", 
                                               y = "Grafting rate (%)") #+ 

#####################
