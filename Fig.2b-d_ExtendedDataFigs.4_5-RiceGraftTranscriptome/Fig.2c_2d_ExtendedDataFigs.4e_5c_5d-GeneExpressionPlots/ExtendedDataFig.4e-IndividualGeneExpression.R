##############################################
##############################################
# Rice graft transcriptome
# R-script by Gregory Reeves (27 MAY 2021)
# University of Cambridge, Department of Plant Sciences,
# Downing Site, Cambridge CB2 3EA, United Kingdom
##############################################
##############################################

library(plyr)
library(tidyverse)
library(stringr)

#New working directory
# THIS WORKING DIRECTORY WILL DEPEND ON YOUR COMPUTER
setwd("~/Documents/CURRENT_PUBLICATIONS/Reeves2020_MonocotGrafting/")


###############################################
###############################################
# ALL TRANSCRIPTOME ANALYSIS HAS BEEN SIMPLIFIED TO USE ONLY ONE INPUT FILE: "fpkm_genename.tsv" 
###############################################
###############################################


### READ IN THE DATA FROM ALL GENE EXPRESSION DATA IN "fpkm_genename.tsv" SAME DATA AS IN "Extended Data Set 2" for R users
fpkm_genename_WIDE <- read.delim("fpkm_genename.tsv", sep="\t")
# head(fpkm_genename_WIDE) # View the data

### CONVERTS DATA INTO LONG FORMAT
fpkm_genename_LONG <- fpkm_genename_WIDE %>% gather(Sample, FPKM, G1D1wc1:G1D5sg6)
# head(fpkm_genename_LONG) # View the data

# This makes a new column called "group" from the Sample name, but removes the final rep number
fpkm_genename_LONG$Group = substr(fpkm_genename_LONG$Sample,1,nchar(fpkm_genename_LONG$Sample)-1)

# This makes a new column called "DAG" Days after grafting from the Sample name, but only keeps day number
fpkm_genename_LONG$DAG = substr(fpkm_genename_LONG$Sample,start = 4, stop = 4)


# This makes a new column called "Treatment" (i.e., ng, wc, sg) from the Sample name
fpkm_genename_LONG$Treatment = substr(fpkm_genename_LONG$Sample,start = 5, stop = 6)

#####################
# We need to remove the G15 (15 days after germination non-grafted samples from the data) from DAG
RemovedGroup <- c("G15ngc")
fpkm_genename_LONG <- filter(fpkm_genename_LONG, !(Group %in% RemovedGroup))

### This converts DAG from a numerical variable into a factor variable
fpkm_genename_LONG$DAG <- as.numeric(fpkm_genename_LONG$DAG)
#####################



###############################################
###############################################
#START HERE FOR ACTUALLY PLOTTING INDIVIDUAL GENES
###############################################
###############################################



#####################
library("ggpubr")

# Filter the data for each gene specifically
# This keeps only one gene 

####CHANGE THE GENE ID HERE TO WHATEVER YOU WANT:
GeneOI <- c("Os05g0194500") 
####

### THIS REMOVES THE WOUNDED CONTROL DATA FOR THE PLOT, optional
# RemoveThisTreatment <- c("wc")

GeneOIfiltered <- filter(fpkm_genename_LONG, (gene_id %in% GeneOI)) # %>% ### Uncomment if you want to filter the wounded control
  # filter(Treatment !=RemoveThisTreatment)  ### Uncomment if you want to filter the wounded control

### This computes a Two-way ANOVA for the gene of interest, and yields the p-value for the Graft type only (not the timepoint)
### Requires library(pylr)
GeneOIpval <- signif(summary(aov(GeneOIfiltered$FPKM ~ GeneOIfiltered$Treatment + GeneOIfiltered$DAG))[[1]][1,5],
                     digits = 4)


### This converts the p-value from above into an asterisk for the plot below. 
GeneOIpvalSignif <- symnum(GeneOIpval, 
                           corr = FALSE, 
                           na = FALSE, 
                           cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 1), 
                           symbols = c("****", "***", "**", "*", ".", "ns"))
GeneOIpvalSignif

#####################
### THIS GENERATES THE PLOTS USED IN THE PUBLICATION
#####################
GeneOILinePlot <- ggline(GeneOIfiltered, x = "DAG", 
                         y = "FPKM", 
                         color = "Treatment", 
                         linetype="Treatment",
                         add = c("mean_se"),   #, "jitter"), ### JITTER CAN BE ADDED TO SEE THE DATA POINTS ON THE PLOT
                         palette = c("#ff754ff8", "#0088aaff", "#7C4700")) # THREE COLORS ARE SPECIFIED BUT ONLY THE FIRST TWO WILL BE USED, For all three treatments


GeneOILinePlot + theme(axis.line = element_line(linetype = "solid"), 
                       axis.ticks = element_line(), 
                       axis.ticks.length = unit(-0.25 , "cm"), ### NEGATIVE length = inside ticks
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
  labs(x = "Days after grafting", y = "FPKM", 
       colour = paste(GeneOI, GeneOIpvalSignif, sep=" ")) +
  annotate("text", label = paste(GeneOIfiltered$gene_name, GeneOI, GeneOIpvalSignif, sep=" "), x = 2.5,
           y=max(GeneOIfiltered$FPKM))




