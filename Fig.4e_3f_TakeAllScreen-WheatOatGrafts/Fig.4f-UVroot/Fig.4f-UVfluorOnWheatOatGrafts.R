##############################################
##############################################
# R-script by Gregory Reeves (27 MAY 2021)
# University of Cambridge, Department of Plant Sciences,
# Downing Site, Cambridge CB2 3EA, United Kingdom
# E-mail me at drgregreeves@gmail.com or gr360@cam.ac.uk
##############################################
##############################################


library(tidyverse)
library(gridExtra)
library(cowplot)
library(ggpubr)
library(multcompView)
library(lsmeans)
library(dplyr)
library(multcompView)
library(agricolae)
library(car) 



### You need to specify the working directory. 
### You can always use the command "file.choose()", then copy the working directory with the file "ImageJ_AnalysisBlueSignalRootArea.csv"
setwd("~/")


rawData <- read.csv("ImageJ_AnalysisBlueSignalRootArea.csv")


#### ONE-WAY ANOVA ANALYSIS


### Equal variance testing (requires car package). If the p-value is greater than 0.05; the variances are homogeneous.
leveneTest(rawData$RelativeInDen , rawData$Plant, center = median)

### Equal variance testing for non-normally distributed data. 
# fligner.test(RelativeInDen ~ Plant, rawData)


### ANOVA TIME
### 

res.aov3 <- aov(RelativeInDen ~ Plant, data = rawData)
summary(res.aov3)

### Summary Stats--Compute mean and SD by groups
group_by(rawData, Plant) %>%
  summarise(
    count = n(),
    mean = mean(RelativeInDen, na.rm = TRUE),
    sd = sd(RelativeInDen, na.rm = TRUE)
  )

# Plot the ANOVA residuals, look for outliers which can affect homogeneity of variance and normality
# It can be useful to remove outliers to meet test assumptions.
plot(res.aov3, 1)

# QQ plot of the ANOVA Resdiuals. 
# The normal probability plot of the residuals should approximately follow a straight line.
plot(res.aov3, 2)

OrderList <- c("oat_selfgraft","oatroot_wheatscion","wheat_selfgraft")


AccessionRename <- c("oat_selfgraft"="Oat-Scion & Oat-Root ",
                     "wheat_selfgraft"="Wheat-Scion & Wheat-Root",
                     "oatroot_wheatscion"="Wheat-Scion & Oat-Root")

### THIS SCRIPT WILL COMPUTE THE ANOVA AGAIN IN ORDER TO EXTRACT TUKEY LETTERS AND DISPLAY THEM ON THE PLOT
#######################################################################

ANOVA_InDen <- aov(RelativeInDen ~ Plant, data=rawData)
TUKEY_InDen <- TukeyHSD(ANOVA_InDen, ordered = FALSE, conf.level = 0.95)


# I need to group the treatments that are not different each other together.
generate_label_df <- function(TUKEY_InDen, variable){
  
  # Extract labels and factor levels from Tukey post-hoc 
  Tukey.levels <- TUKEY_InDen[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  
  #I need to put the labels in the same order as in the boxplot :
  Tukey.labels$Plant=rownames(Tukey.labels)
  Tukey.labels=Tukey.labels[order(Tukey.labels$Plant) , ]
  return(Tukey.labels)
}

# Apply the function on my dataset
LABELSInDen=generate_label_df(TUKEY_InDen , "Plant")

#Reorder the TUKEY letters to match the order of the accessions in the plots 
LABELSInDennew <-LABELSInDen %>%
  slice(match(OrderList, Plant))

#Get the max values from each accession so that the letters know where to go

meansInDen <- aggregate(RelativeInDen ~  Plant, rawData, max)

meansPlantnew <-meansInDen %>%
  slice(match(OrderList, Plant))



#Generate the plot itself 

AnInDen <- ggboxplot(rawData, 
                     x="Plant", 
                     y="RelativeInDen", 
                     #y="IntDen_blue",
                     order=OrderList, 
                     fill= "Plant", 
                     color = "black", 
                     weight= 3, 
                     add = "jitter",  
                     #size = 1,
                     add.params = list(color = "Plant", size=2))+
  scale_y_continuous(name=expression("Normalized blue fluorescence", expand=c(0.001,0.01))) +
  theme_classic() + theme(legend.title=element_blank()) +
  scale_fill_manual(breaks = c("oat_selfgraft", "wheat_selfgraft", "oatroot_wheatscion"), 
                    values=c("#2ca089ff", "#502d16ff", "#677821ff")) +
  scale_color_manual(breaks = c("oat_selfgraft", "wheat_selfgraft", "oatroot_wheatscion"), 
                     values=c("#37c8abff", "#504416ff", "#88aa00ff")) +
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
  labs(x = "Graft Type") + 
  stat_compare_means(label = "p.signif", method = "anova", label.y =0.7, label.x=2, size = 8)+
  scale_x_discrete(labels=AccessionRename)
#  Add p-value
AnInDen + annotate("text", x = 1:3, 
                   y = meansPlantnew$RelativeInDen,
                   label = LABELSInDennew$Letters,
                   vjust=-1.5,
                   fontface =3,
                   size = 6.5)




