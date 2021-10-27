
###################################################
# Monocot grafting
# R-script by Gregory Reeves (11 July 2021)
# gr360@cam.ac.uk
# University of Cambridge, Department of Plant Sciences
###################################################

library(tidyverse)
library(dplyr) # Used for the ANOVA analysis 
library(car) # Used for Levene's test
library(multcompView)
library(ggpubr)
library(ggplot2)
library(viridis)

### READ IN THE DATA
### You need to specify the working directory. 
### You can always use the command "file.choose()", then copy the working directory with the file "RiceTimecourseTranscriptome-21Oct2020.csv"
#setwd("~/")

file.choose() ### You need to choose the file called "SEMcellLength.csv"

##############################
### Cell_Length Plots
##############################

#Grafting data input
SEMCellLength <- read.csv("SEMcellLength.csv", sep=",",header=TRUE)

SEMCellLength <- SEMCellLength %>% mutate(ComboDAG = paste(Combo, DAG, sep="_"))

# Rename Columne in AttachmentData
colnames(SEMCellLength) <-c("DAG",	
                            "Half",	
                            "Type",	
                            "Combo",
                            "CellLength",
                            "PlantRep",
                            "ComboDAG")

# Summarise data
SEMCellLengthSummary <- summarise(group_by(SEMCellLength, 
                                           DAG, Half, Type, Combo, PlantRep, ComboDAG),
                                  meanCellLength=mean(CellLength),
                                  sdCellLength=sd(CellLength),
                                  seCellLength=sd(CellLength)/sqrt(n()),
                                  n=n())


####################
# ANOVA ANALYSIS
###################


###Remove nongrafted controls, if desired 
Remove_ngc <- c("Nongrafted")
SEMCellLengthSummaryGraftOnly <- filter(SEMCellLengthSummary, !(Type %in% Remove_ngc))

#SEMCellLengthSummaryGraftOnly <- SEMCellLengthSummary ###Uncomment/Use if not filtering

# Check the structure
str(SEMCellLengthSummaryGraftOnly)

# Convert DAG as a factor and recode the levels
# as "1DAG", "3DAG", "5DAG", "7DAG"
SEMCellLengthSummaryGraftOnly$DAG <- factor(SEMCellLengthSummaryGraftOnly$DAG, 
                                            levels = c(1, 3, 5, 7),
                                            labels = c("1DAG", "3DAG", "5DAG", "7DAG"))
head(SEMCellLengthSummaryGraftOnly)


### Generate frequency tables, We want to know if Days After Grafting, and which graft type affect CellLength
table(SEMCellLengthSummaryGraftOnly$DAG, 
      SEMCellLengthSummaryGraftOnly$Combo)

### Result: We have a 4X4 design, balanced with 5 observations in each cell

# Box plot with multiple groups
# +++++++++++++++++++++
# Plot Cell length ("meanCellLength") by days after grafting("DAG")
# Color box plot by a second group, root or scion cells: "Combo"

library("ggpubr")
ggboxplot(SEMCellLengthSummaryGraftOnly, 
          x = "DAG", 
          y = "meanCellLength", 
          color = "Combo",
          palette = c("#00AFBB", "#E7B800"))

### library(car)
### Check the homogeneity of variances
#leveneTest(meanCellLength ~ DAG*Combo, data = SEMCellLengthSummaryGraftOnly)

CellLengthANOVA <- aov(meanCellLength ~ Combo + DAG + Combo:DAG, data = SEMCellLengthSummaryGraftOnly)
summary(CellLengthANOVA) ### THIS GENERATES THE ANOVA TABLE AND THE P-VALES IN Fig.1k

#### From the ANOVA table we can conclude that both Combo and DAG are statistically significant, as well as their interaction. 
#### Note the above fitted model is called additive model. It makes an assumption that the two factor variables are independent. 
#### If you think that these two variables might interact to create an synergistic effect, replace the plus symbol (+) by an asterisk (*), or by adding "+ Combo:DAG". 
#### Here we can see that there is interaction between DAG and Combo. Note that, in the situation where the interaction is not significant you should use the additive model.

TukeyHSD(CellLengthANOVA, which = "Combo:DAG")

plot(TukeyHSD(CellLengthANOVA, which = "Combo:DAG"))

Output <- TukeyHSD(CellLengthANOVA, which = "Combo:DAG")
PT4 <- Output$Combo[,'p adj']
multcompLetters(PT4) ### THIS YIELDS THE SIGNIFICANT GROUPINGS ACCORDING TO THE TUKEY TEST SHOWN IN Fig.1k


### Check the homogeneity of variances
# 1. Homogeneity of variances
plot(CellLengthANOVA, 1)

# 2. Normality
plot(CellLengthANOVA, 2)
### As all the points fall approximately along this reference line, we can assume normality.

# Extract the residuals
aov_residuals <- residuals(object = CellLengthANOVA)
# Run Shapiro-Wilk test
shapiro.test(x = aov_residuals)
### The conclusion above, is supported by the Shapiro-Wilk test on the ANOVA residuals (W = 0.98, p = 0.74) which finds slight indication that normality is violated (near to  p = 0.05).


####################
# Box plot (with scatter overlay)
###################

###Remove nongrafted controls, if desired 
Remove_ngc <- c("Nongrafted")
SEMCellLengthGraftOnly <- filter(SEMCellLength, !(Type %in% Remove_ngc))


ggplot(data = SEMCellLengthSummaryGraftOnly,           ## for n = 5, use: data = SEMCellLengthSummaryGraftOnly, and y = meanCellLength,
       aes(x = ComboDAG, y = meanCellLength, fill = Combo)) + 
  #geom_jitter() +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) +
  geom_violin(alpha = 0.5, width = 1, trim = FALSE,  adjust = 1) +
  geom_boxplot(outlier.color = "black", alpha = 0.6, width = 0.15) +
  #scale_x_continuous(breaks=seq(1,7,by=2)) +
  scale_fill_manual(values = c("#95731d", "#aeba0e")) +
  #facet_grid(. ~DAG) +
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
                                   angle = 90,
                                   margin=margin(10,5,10,5,"pt")), 
        axis.text.y = element_text(size = 15,
                                   #angle = 45,
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
  labs(x = "Days after grafting", y = "Average cell length (µm)")

