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
#######################

library(ggplot2)
library(tidyverse)
library(viridis)
#library(Hmisc) #Allows minor.ticks function in ggplot


#####################
# Other Order Grafts
#####################



### You need to specify the working directory. 
### You can always use the command "file.choose()", then copy the working directory with the file "OtherOrdersGrafts.csv"
setwd("~/")




#Grafting data input
OtherOrders <- read.csv("OtherOrdersGrafts.csv", sep=",",header=TRUE)


OtherOrders_summarised <- OtherOrders %>%
  group_by(Order, 
           Species) %>%
  summarise(totalFused=sum(Fused),
            totalAttempts=sum(Attempts),
            FusionRate=100*(sum(Fused)/sum(Attempts)))

OtherOrders_summarised <- OtherOrders_summarised[-c(1), ] #Removes row 1 which appears to be messed up

#####################
#We need to filter by each ORDER to make each plot
OtherOrdersFILTERED <-  filter(OtherOrders_summarised, (Order == 'Alismatales'))

###Poales
#OrderLevels <- c('Triticum aestivum', 'Triticum durum', 'Secale cereale','Avena sativa', 
#                 'Phyllostachys edulis', 'Oryza sativa', 'Pennisetum glaucum', 'Zea mays', 
#                 'Sorghum bicolor', 'Puya alpestris', 'Puya raimondii', 'Ananas comosus') 

###Zingiberales
#OrderLevels <- c('Musa acuminata', 'Musa balbisiana', 'Amomum sublatum', 'Canna indica', 'Costus laevis') 

###Commelinales
#OrderLevels <- c('Pollia japonica', 'Pollia thysiflora', 'Commelina comminis') 

###Areceales
#OrderLevels <- c('Phoenix dactylifera', 'Elaeis guineensis', 'Elaeis oleifera') 

###Asparagales
#OrderLevels <- c('Allium cepa', 'Allium schoenoprasum', 'Allium ampeloprasum', 'Beaucarnea recurvata', 'Phormium tenax', 'Yucca brevifolia', 'Agave tequilana') 


###Dioscoreales
#OrderLevels <- c('Tacca leontopetaloides', 'Dioscorea elephantipes') 

###Liliales
#OrderLevels <- c('Bomarea multiflora', 'Gloriosa superba') 

###Alismatales
#OrderLevels <- c('Arisaema tortuosum') 

###Acorales
#OrderLevels <- c('Acorus calamus') 


###ALL ORDER SPECIES ORDERED
OrderLevels <- c('Triticum aestivum', 'Triticum durum', 'Secale cereale','Avena sativa', 
                 'Phyllostachys edulis', 'Oryza sativa', 'Pennisetum glaucum', 'Zea mays', 
                 'Sorghum bicolor', 'Puya alpestris', 'Puya raimondii', 'Ananas comosus',
                 'Musa acuminata', 'Musa balbisiana', 'Amomum sublatum', 'Canna indica', 'Costus laevis',
                 'Pollia japonica', 'Pollia thysiflora', 'Commelina comminis',
                 'Phoenix dactylifera', 'Elaeis guineensis', 'Elaeis oleifera',
                 'Allium cepa', 'Allium schoenoprasum', 'Allium ampeloprasum', 'Beaucarnea recurvata', 'Phormium tenax', 'Yucca brevifolia', 'Agave tequilana',
                 'Tacca leontopetaloides', 'Dioscorea elephantipes',
                 'Bomarea multiflora', 'Gloriosa superba',
                 'Arisaema tortuosum',
                 'Acorus calamus')




#####################

### Other Orders PLOT


ggplot(data = OtherOrdersFILTERED, 
       aes(x=factor(Species, levels = OrderLevels),
           y=FusionRate,
           fill = Species)) + 
  
  #geom_col(position = position_dodge2(width = 0.9, preserve = "single")) +
  geom_bar(stat="identity", position="dodge") +
  geom_text(aes(label = paste("n =", totalAttempts, sep=" "), 
                vjust = -1),
            color = "gray30",
            size = 4) + 
  scale_y_continuous(breaks=seq(0,80,by=5),
                     limits=c(0,80)) +
  #scale_fill_viridis_d() +
  scale_fill_grey() +
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
        legend.direction = "horizontal") +labs(x = "Species", 
                                               y = "Grafting rate (%)") #+ 

#####################





#####################

### Other Orders PLOT


ggplot(data = OtherOrders_summarised, 
       aes(#x=Species,
         x=factor(Species, levels = OrderLevels),
         y=FusionRate,
         fill = Species)) + 
  geom_col(position = position_dodge2(width = 0.9, preserve = "single")) +
  #geom_bar(stat="identity", position="dodge") +
  #geom_text(aes(label = paste("n =", totalAttempts, sep=" "), vjust = -1), color = "gray30", size = 4) + 
  geom_text(aes(label = paste(totalAttempts), vjust = -1), color = "gray30", size = 4) + 
  scale_y_continuous(breaks=seq(0,80,by=5),
                     limits=c(0,80)) +
  facet_grid(~ Order, scales = "free_x", space = "free_x") + # Comment this out to intercalate the bars
  #scale_fill_viridis_d() +
  scale_fill_grey() +
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
                                   angle = 90,
                                   face = "italic"), 
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
        legend.direction = "horizontal") +labs(x = "Species", 
                                               y = "Grafting rate (%)") #+ 

#####################


