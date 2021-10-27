##############################################
##############################################
# R-script by Gregory Reeves 30 April 2020
# University of Cambridge, Department of Plant Sciences,
# Downing Site, Cambridge CB2 3EA, United Kingdom
# E-mail me at drgregreeves@gmail.com 
##############################################
##############################################

library(tidyverse)
library(gridExtra)
library(cowplot)
library(ggpubr)
library(dbplyr)
library(rstatix) ## Required for the repeated measures two-way ANOVA
# library(ggThemeAssist) # Usage ggThemeAssistGadget(plotname)



### You need to specify the working directory. 
### You can always use the command "file.choose()", then copy the working directory with the file "ImageJ_AnalysisBlueSignalRootArea.csv"
#setwd("~/")


StrigolactoneRiceMutantGrafting <- file.choose() ### Pick the file called "StrigolactoneRiceMutantGrafting_allmerged.csv"


StrigoData <- read.csv(StrigolactoneRiceMutantGrafting)



colnames(StrigoData) <-c("id", "DaysOld", "Graft", "Plant", "Tiller","Ligule", "Leaf")

StrigoDataSummary <- summarise(group_by(StrigoData, DaysOld, Graft),
                               meanTil=mean(Tiller),
                               meanLigule=mean(Ligule),
                               meanLeaf=mean(Leaf),
                               sdTil=sd(Tiller),
                               sdLigule=sd(Ligule),
                               sdLeaf=sd(Leaf),
                               seTil=sd(Tiller)/sqrt(n()),
                               seLigule=sd(Ligule)/sqrt(n()),
                               seLeaf=sd(Leaf)/sqrt(n()),
                               n=n())


###########
#### Fig. 4b -- code for the line plot
###########

# With standard deviation
TillerPlot <- ggplot(data=StrigoDataSummary, mapping=aes(DaysOld,meanTil,group=Graft, color=Graft)) +
  geom_line(size=2) + 
  geom_errorbar(aes(DaysOld, ymin=meanTil-sdTil, ymax=meanTil+sdTil, width=2)) +
  #geom_jitter(data=StrigoData, aes(DaysOld, Tiller, shape=Graft, color=Graft, fill=Graft), inherit.aes = FALSE)+
  geom_point(aes(shape=Graft, color=Graft, fill=Graft), size=4)+ 
  scale_shape_manual(values=c(21,22,23))+ 
  scale_fill_manual(breaks = c("Mu-Mu", "Wt-Mu", "Wt-Wt"), values=c("#1c344eff", "#92c46dff", "#297d7dff"))+
  labs(x = expression("Days after grafting"), y = expression("No. tillers")) +   
  scale_color_manual(breaks = c("Mu-Mu", "Wt-Mu", "Wt-Wt"), values=c("#485d70ff","#a7d189ff", "#539798ff")) +
  scale_x_continuous(breaks=seq(0,60,5))+
  scale_y_continuous(breaks=seq(0, 50, 5))+
  theme(legend.position="none",
        axis.line.x = element_line(colour = "#000000", size = 2, linetype = "solid"),
        axis.line.y = element_line(colour = "#000000", size = 2, linetype = "solid"),
        axis.ticks = element_line(colour = "#000000", size = 2, linetype = "solid"),
        axis.text.x = element_text(face="bold", size=30, colour = "#333333"),
        axis.text.y = element_text(face="bold", size=30, colour = "#333333"),
        axis.title = element_text(size=30)) 
TillerPlot + theme(panel.grid.major = element_line(colour = NA), 
                   panel.grid.minor = element_line(colour = NA), 
                   panel.background = element_rect(fill = NA))



#########
### Extended Data Fig. 9b -- code for the line plots
#########

# With standard deviation
LigulePlot <- ggplot(data=StrigoDataSummary, mapping=aes(DaysOld,meanLigule,group=Graft, color=Graft)) +
  geom_line(size=2, position=position_dodge(width=0.2)) + 
  geom_errorbar(aes(DaysOld, ymin=meanLigule-sdLigule, ymax=meanLigule+sdLigule, width=2)) +
  #geom_jitter(data=StrigoData, aes(DaysOld, Ligule, shape=Graft, color=Graft, fill=Graft), inherit.aes = FALSE)+
  geom_point(aes(shape=Graft, color=Graft, fill=Graft), size=4)+ 
  scale_shape_manual(values=c(21,22,23))+ 
  scale_fill_manual(breaks = c("Mu-Mu", "Wt-Mu", "Wt-Wt"), values=c("#1c344eff", "#92c46dff", "#297d7dff"))+
  labs(x = expression("Days after grafting"), y = expression("Tallest ligule (cm)")) +   
  scale_color_manual(breaks = c("Mu-Mu", "Wt-Mu", "Wt-Wt"), values=c("#485d70ff","#a7d189ff", "#539798ff")) +
  scale_x_continuous(breaks=seq(0,60,5))+
  scale_y_continuous(breaks=seq(0, 60, 5))+
  theme(legend.position="none",
        axis.line.x = element_line(colour = "#000000", size = 2, linetype = "solid"),
        axis.line.y = element_line(colour = "#000000", size = 2, linetype = "solid"),
        axis.ticks = element_line(colour = "#000000", size = 2, linetype = "solid"),
        axis.text.x = element_text(face="bold", size=30, colour = "#333333"),
        axis.text.y = element_text(face="bold", size=30, colour = "#333333"),
        axis.title = element_text(size=30)) 
LigulePlot + theme(panel.grid.major = element_line(colour = NA), 
                   panel.grid.minor = element_line(colour = NA), 
                   panel.background = element_rect(fill = NA))


# With standard deviation
LeafPlot <- ggplot(data=StrigoDataSummary, mapping=aes(DaysOld,meanLeaf,group=Graft, color=Graft)) +
  geom_line(size=2, position=position_dodge(width=0.2)) + 
  geom_errorbar(aes(DaysOld, ymin=meanLeaf-sdLeaf, ymax=meanLeaf+sdLeaf, width=2)) +
  #geom_jitter(data=StrigoData, aes(DaysOld, Leaf, shape=Graft, color=Graft, fill=Graft), inherit.aes = FALSE)+
  geom_point(aes(shape=Graft, color=Graft, fill=Graft), size=4)+ 
  scale_shape_manual(values=c(21,22,23))+ 
  scale_fill_manual(breaks = c("Mu-Mu", "Wt-Mu", "Wt-Wt"), values=c("#1c344eff", "#92c46dff", "#297d7dff"))+
  labs(x = expression("Days after grafting"), y = expression("Longest leaf (cm)")) +   
  scale_color_manual(breaks = c("Mu-Mu", "Wt-Mu", "Wt-Wt"), values=c("#485d70ff","#a7d189ff", "#539798ff")) +
  scale_x_continuous(breaks=seq(0,60,5))+
  scale_y_continuous(breaks=seq(0, 110, 5))+
  theme(legend.position="none",
        axis.line.x = element_line(colour = "#000000", size = 2, linetype = "solid"),
        axis.line.y = element_line(colour = "#000000", size = 2, linetype = "solid"),
        axis.ticks = element_line(colour = "#000000", size = 2, linetype = "solid"),
        axis.text.x = element_text(face="bold", size=30, colour = "#333333"),
        axis.text.y = element_text(face="bold", size=30, colour = "#333333"),
        axis.title = element_text(size=30)) 
LeafPlot  + theme(panel.grid.major = element_line(colour = NA), 
                  panel.grid.minor = element_line(colour = NA), 
                  panel.background = element_rect(fill = NA))


####################################
### Use this if you want all three plots together
####################################

StrigoFig <- plot_grid( TillerPlot, LeafPlot, LigulePlot,
                        ncol = 3, nrow = 1, align = "hv", labels = "AUTO")
StrigoFig




data <- StrigoData

### Plot the data
### requires "ggpubr"

#########################
### Same plot as Fig. 4b
#########################
TillerLinePlot <- ggline(data, x = "DaysOld", y = "Tiller", color = "Graft",
                         add = c("mean_sd", "jitter"),
                         palette = c("#1c344eff", "#92c46dff", "#297d7dff"))
TillerLinePlot


TillerLinePlot + theme(axis.line = element_line(size = 1, 
                                                linetype = "solid"), axis.ticks = element_line(size = 1), 
                       panel.grid.major = element_line(colour = NA), 
                       panel.grid.minor = element_line(colour = NA), 
                       axis.title = element_text(size = 18, 
                                                 face = "bold"), axis.text = element_text(size = 15, 
                                                                                          face = "bold", colour = "gray30"), 
                       axis.text.x = element_text(size = 15), 
                       axis.text.y = element_text(size = 15), 
                       legend.text = element_text(size = 9), 
                       legend.title = element_text(size = 15, 
                                                   face = "bold", colour = "gray30"), 
                       panel.background = element_rect(fill = NA), 
                       legend.key = element_rect(fill = NA), 
                       legend.background = element_rect(fill = NA), 
                       legend.direction = "horizontal") +labs(x = "Days after grafting", y = "Tiller (no.)", 
                                                              colour = "Graft type ")



################
################
##### two-way ANOVA repeated measures ANOVA on DATA
#--------------------------------


#########################
# For TILLER
#########################
data <- StrigoData

set.seed(1234)
dplyr::sample_n(data, 10) # See the data in a table format

str(data) #Seeing what R considers each variable type to be

### This converts DaysOld from a numerical variable into a factor variable
data$DaysOld <- factor(data$DaysOld, 
                       levels = c(8, 14, 36, 43, 55),
                       labels = c("8", "14", "36", "43", "55"))
head(data)


data <- StrigoData
data$id <- as.factor(data$id)
# data$DaysOld <- as.factor(data$DaysOld)
data$Plant <- as.factor(data$Plant)
# data <- data %>% select(id, DaysOld, Graft, Tiller)
data <- data %>% select(Plant, DaysOld, Graft, Tiller)

# data <- data %>% filter(Graft %in% c("Mu-Mu", "Wt-Wt"))
# data <- data %>% filter(Tiller > 1)
data <- data %>% filter(DaysOld > 14)
data$DaysOld <- as.factor(data$DaysOld)

# Inspect some random rows of the data by groups
set.seed(123)
data %>% sample_n_by(Graft, DaysOld, size = 1)


# Summary Stats
data %>%
  group_by(Graft, DaysOld) %>%
  get_summary_stats(Tiller, type = "mean_sd")

bxp <- ggboxplot(
  data, x = "DaysOld", y = "Tiller",
  color = "Graft", palette = "jco"
)
bxp

# Outliers can be easily identified using this:
data %>%
  group_by(Graft, DaysOld) %>%
  identify_outliers(Tiller)

#The normality assumption can be checked
#by computing Shapiro-Wilk test for each time point.
#If the data is normally distributed, 
#the p-value should be greater than 0.05.
data %>%
  group_by(Graft, DaysOld) %>%
  shapiro_test(Tiller)

# Or by QQ plot

qqplotTiller <- ggqqplot(data, "Tiller", ggtheme = theme_bw()) +
  facet_grid(DaysOld ~ Graft, labeller = "label_both")
qqplotTiller 

#### Equal variance testing

# data %>% levene_test(Tiller ~ Graft*DaysOld)


# Run the actual ANOVA

res.aov <- anova_test(
  data = data, dv = Tiller, wid = Plant,
  within = c(Graft, DaysOld)
)
get_anova_table(res.aov) ### THE TWO-WAY REPEATED MEASURES ANOVA TABLE



# Effect of treatment at each time point
one.way <- data %>%
  group_by(DaysOld) %>%
  anova_test(dv = Tiller, wid = Plant, within = Graft) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")
one.way

pwc2 <- data %>%
  group_by(DaysOld) %>%
  pairwise_t_test(
    Tiller ~ Graft, paired = TRUE,
    p.adjust.method = "bonferroni"
  )
pwc2



#########################
# For LIGULE
#########################
data <- StrigoData

set.seed(1234)
dplyr::sample_n(data, 10) # See the data in a table format

str(data) #Seeing what R considers each variable type to be

### This converts DaysOld from a numerical variable into a factor variable
data$DaysOld <- factor(data$DaysOld, 
                       levels = c(8, 14, 36, 43, 55),
                       labels = c("8", "14", "36", "43", "55"))
head(data)


data <- StrigoData
data$id <- as.factor(data$id)
# data$DaysOld <- as.factor(data$DaysOld)
data$Plant <- as.factor(data$Plant)
# data <- data %>% select(id, DaysOld, Graft, Tiller)
data <- data %>% select(Plant, DaysOld, Graft, Ligule)

# data <- data %>% filter(Graft %in% c("Mu-Mu", "Wt-Wt"))
# data <- data %>% filter(Tiller > 1)
# data <- data %>% filter(DaysOld > 14)
data$DaysOld <- as.factor(data$DaysOld)

# Inspect some random rows of the data by groups
set.seed(123)
data %>% sample_n_by(Graft, DaysOld, size = 1)


# Summary Stats
data %>%
  group_by(Graft, DaysOld) %>%
  get_summary_stats(Ligule, type = "mean_sd")

bxp <- ggboxplot(
  data, x = "DaysOld", y = "Ligule",
  color = "Graft", palette = "jco"
)
bxp

# Outliers can be easily identified using this:
data %>%
  group_by(Graft, DaysOld) %>%
  identify_outliers(Ligule)

#The normality assumption can be checked
#by computing Shapiro-Wilk test for each time point.
#If the data is normally distributed, 
#the p-value should be greater than 0.05.
data %>%
  group_by(Graft, DaysOld) %>%
  shapiro_test(Ligule)

# Or by QQ plot

qqplotLigule <- ggqqplot(data, "Ligule", ggtheme = theme_bw()) +
  facet_grid(DaysOld ~ Graft, labeller = "label_both")
qqplotLigule
#### Equal variance testing

# data %>% levene_test(Ligule ~ Graft*DaysOld)


res.aov <- anova_test(
  data = data, dv = Ligule, wid = Plant,
  within = c(Graft, DaysOld)
)
get_anova_table(res.aov) ### THE TWO-WAY REPEATED MEASURES ANOVA TABLE


# Effect of treatment at each time point
one.way <- data %>%
  group_by(DaysOld) %>%
  anova_test(dv = Ligule, wid = Plant, within = Graft) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")
one.way




#########################
# For LEAF
#########################
data <- StrigoData

set.seed(1234)
dplyr::sample_n(data, 10) # See the data in a table format

str(data) #Seeing what R considers each variable type to be

### This converts DaysOld from a numerica variable into a factor variable
data$DaysOld <- factor(data$DaysOld, 
                       levels = c(8, 14, 36, 43, 55),
                       labels = c("8", "14", "36", "43", "55"))
head(data)


data <- StrigoData
data$id <- as.factor(data$id)
# data$DaysOld <- as.factor(data$DaysOld)
data$Plant <- as.factor(data$Plant)
# data <- data %>% select(id, DaysOld, Graft, Tiller)
data <- data %>% select(Plant, DaysOld, Graft, Leaf)

# data <- data %>% filter(Graft %in% c("Mu-Mu", "Wt-Wt"))
# data <- data %>% filter(Tiller > 1)
# data <- data %>% filter(DaysOld > 14)
data$DaysOld <- as.factor(data$DaysOld)

# Inspect some random rows of the data by groups
set.seed(123)
data %>% sample_n_by(Graft, DaysOld, size = 1)


# Summary Stats
data %>%
  group_by(Graft, DaysOld) %>%
  get_summary_stats(Leaf, type = "mean_sd")

bxp <- ggboxplot(
  data, x = "DaysOld", y = "Leaf",
  color = "Graft", palette = "jco"
)
bxp

# Outliers can be easily identified using this:
data %>%
  group_by(Graft, DaysOld) %>%
  identify_outliers(Leaf)

#The normality assumption can be checked
#by computing Shapiro-Wilk test for each time point.
#If the data is normally distributed, 
#the p-value should be greater than 0.05.
data %>%
  group_by(Graft, DaysOld) %>%
  shapiro_test(Leaf)

# Or by QQ plot

qqplotLeaf <- ggqqplot(data, "Leaf", ggtheme = theme_bw()) +
  facet_grid(DaysOld ~ Graft, labeller = "label_both")
qqplotLeaf


#### Equal variance testing

#data %>% levene_test(Leaf ~ Graft*DaysOld)


res.aov <- anova_test(
  data = data, dv = Leaf, wid = Plant,
  within = c(Graft, DaysOld)
)
get_anova_table(res.aov)  ### THE TWO-WAY REPEATED MEASURES ANOVA TABLE


# Effect of treatment at each time point
#Note: that only Bonferroni and Sidak (besides LSD) are 
#done for repeated measure factor
one.wayLeaf <- data %>%
  group_by(DaysOld) %>%
  anova_test(dv = Leaf, wid = Plant, within = Graft) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")
one.wayLeaf

pwc2 <- data %>%
  group_by(DaysOld) %>%
  pairwise_t_test(
    Leaf ~ Graft, paired = TRUE,
    p.adjust.method = "bonferroni"
  )
pwc2





####################################
# QQ Plot figure merged (requires cowplot package)

qqFig <- plot_grid(qqplotTiller, qqplotLigule, qqplotLeaf,
                   ncol = 3, nrow = 1, align = "hv", labels = "AUTO")
qqFig


###############################





