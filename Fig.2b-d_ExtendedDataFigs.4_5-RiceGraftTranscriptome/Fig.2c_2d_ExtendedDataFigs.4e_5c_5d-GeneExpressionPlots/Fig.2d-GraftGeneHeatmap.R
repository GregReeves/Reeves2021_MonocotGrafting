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

library(gplots)
library(viridis)


#New working directory
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
#START HERE FOR ACTUALLY PLOTTING HEATMAPS
###############################################
###############################################


#############################
# IF filtering by specific genes
############################

###Remove wounded control (wc) if desired 
Remove_wc <- c("wc")
fpkm_genename_LONG_forZ <- filter(fpkm_genename_LONG, !(Treatment %in% Remove_wc))

### Run this if not filtering out Wounded
#fpkm_genename_LONG_forZ <- fpkm_genename_LONG

#Computes the Z-score for every gene at each timepoint 
fpkm_genename_LONG_Z <- fpkm_genename_LONG_forZ %>% 
  group_by(gene_id, DAG) %>%
  mutate(zscore = (FPKM - mean(FPKM))/sd(FPKM))



###These are the rice gene IDs associated with grafting in Arabidopsis that are correlated in rice, and some novel genes from rice only ###Uncomment with "shift/control/c"
RiceCLOSEGraftingGenesIDs <- c("Os06g0726400",
                               "Os09g0298200",
                               "Os01g0851700",
                               "Os10g0555700",
                               "Os09g0491740",
                               "Os07g0589000",
                               "Os02g0228900",
                               "Os01g0277700",
                               "Os03g0243700",
                               "Os06g0696600",
                               "Os02g0696500",
                               "Os04g0398000",
                               "Os12g0641400",
                               "Os05g0194500",
                               "Os02g0663100",
                               "Os10g0442600",
                               "Os06g0214800",
                               "Os06g0214850",
                               "Os01g0313300",
                               "Os02g0662700",
                               "Os03g0760800",
                               "Os10g0571600",
                               "Os06g0222400",
                               "Os05g0444200",
                               "Os09g0485500",
                               "Os12g0560300",
                               "Os02g0513100",
                               "Os01g0358100",
                               "Os02g0707200",
                               "Os02g0521100",
                               "Os03g0157600",
                               "Os08g0290100",
                               "Os03g0687000",
                               "Os03g0111000",
                               "Os09g0528200",
                               "Os11g0275200",
                               "Os01g0631500",
                               "Os04g0444300",
                               "Os02g0237100",
                               "Os04g0497200",
                               "Os05g0515700",
                               "Os04g0627000",
                               "ENSRNA049471001",
                               "ENSRNA049471381",
                               "ENSRNA049470347")

#Checks for duplicate data 
RiceCLOSEGraftingGenesIDs[duplicated(RiceCLOSEGraftingGenesIDs)] 

#Filters all the expression data for the genes of interest
THELIST <- fpkm_genename_LONG_Z %>%
  filter(gene_id %in% RiceCLOSEGraftingGenesIDs) 




#####  THIS  STEP RUNS A TWO-WAY ANOVA ON ALL GENES SPECIFIED, AND ADDS THE GROUP P-VALUE AS A NEW COLUMN, REQUIRES library(plyr)
ANOVAdForGeneFamilyPlot <- ddply(THELIST,.(gene_id),
                                 transform,
                                 pval= summary(aov(FPKM ~ Treatment + DAG))[[1]][1,5])
# head(ANOVAdForGeneFamilyPlot) # View the data


# THIS REMOVED ANY Na VALUES FROM THE P-VALUE COLUMN WHICH WILL INTERFERE WITH THE HEATMAP
ANOVAdForGeneFamilyPlot <-ANOVAdForGeneFamilyPlot %>% 
  filter(pval != c("NaN")) %>%
  filter(zscore != c("NaN")) 


### WE NEED TO USE SOME FUNCTIONS FROM THE dplyr PACKAGE
detach(package:plyr) # dplyr Does not work with plyr loaded
### NOW WE LOAD dplyr
 library(dplyr)

### THIS COMPUTES THE AVERAGE P-VALUE/FPKM/Z-SCORE/ETC
GeneFamilyfilteredSummary <- ANOVAdForGeneFamilyPlot %>%
  group_by(gene_id, Treatment, DAG) %>%
  summarise(meanFPKM=mean(FPKM), 
            meanZscore=mean(zscore),
            sdFPKM=sd(FPKM), 
            seFPKM=sd(FPKM)/sqrt(3), 
            mean_pval=mean(pval))

#### OPTIONAL: CHANGE THE mean_pval LEVEL TO FILTER ON P-VALUE, HERE WE ARE NOT FILTERING (anything less than 1). 
GeneFamilyfilteredSummary <-GeneFamilyfilteredSummary %>% 
  filter(mean_pval < 1)

###Combines the treatment and timepoint into new column
GeneFamilyfilteredSummary$treatmentDAG <- paste0(GeneFamilyfilteredSummary$Treatment, 
                                                 GeneFamilyfilteredSummary$DAG)

###Removes the original treatment, timepoint, sd, se, and p-value columns
GeneFamilyfilteredSummary <- GeneFamilyfilteredSummary %>% 
  ungroup() %>% 
  select(-DAG, -Treatment, -mean_pval, -sdFPKM, -seFPKM, -meanFPKM)

###Converts the entire dataframe into wide format (to help create the matrix) based on Z score
GeneFamilyfilteredSummary_wide <- GeneFamilyfilteredSummary %>% 
  pivot_wider(names_from = treatmentDAG, 
              values_from = meanZscore)


###This formats the matrix so that it can be used for the heatmap function
df_matrix <- as.matrix(select(GeneFamilyfilteredSummary_wide, -gene_id))
rownames(df_matrix) <- GeneFamilyfilteredSummary_wide$gene_id

###This generates the heatmap itself for the specific search term
heatmap.2(df_matrix,
          Colv = F, 
          #col = brewer.pal(n = 8, name = "Greens"), # Uncomment to use Brewer Colors
          col = viridis_pal(alpha = 1, begin = 1, end = 0, direction = 1, option = "D"),# Uncomment to use Viridis Colors 
          dendrogram = "row", 
          trace = "none", 
          symkey = F, 
          scale = "none", 
          cexRow = 1, 
          key = T, 
          key.title = NA, 
          keysize = 1)


### enjoy your hot heatmap ;)














