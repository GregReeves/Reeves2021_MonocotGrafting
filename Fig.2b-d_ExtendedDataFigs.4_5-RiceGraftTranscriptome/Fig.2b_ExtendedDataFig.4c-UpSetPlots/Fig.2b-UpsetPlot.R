##############################################
##############################################
# Rice graft transcriptome
# R-script by Gregory Reeves (27 MAY 2021)
# University of Cambridge, Department of Plant Sciences,
# Downing Site, Cambridge CB2 3EA, United Kingdom
# E-mail me at drgregreeves@gmail.com or gr360@cam.ac.uk
##############################################
##############################################


library(tidyverse)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(UpSetR) # TO PLOT AS UPSET DIAGRAMS
library(VennDiagram) # TO PLOT AS VENN DIAGRAMS

# Set parent directory -- THIS NEEDS TO SPECIFY THE FOLDER "coExpression_venn"
setwd("~/")


#### #######
#### THIS IS A CUSTOM FUNCTION CALLED "readDataFiles" THAT EXTRACTS THE CORRECT EXPRESSION DATA FROM THE "coExpression_venn" FOLDER
#### #######

readDataFiles <- function(folder) {
  fileList <- list.files(pattern="*.tsv", path=folder, full.names=TRUE)
  output <- lapply(fileList, read_tsv)
  names(output) <- gsub(pattern=paste(folder, "/", sep=""), replacement="", fileList)
  output <- lapply(output, function(df) mutate_at(df, .vars = c("gene_chr"), as.character))
  return(output)
}




######## Fig. 2b
# LETS EXTRACT grafted VS non-grafted UP-regulated genes only
######## 
G1D1sgvsG1D1ngc_G1D1wcvsG1D1ngc <- readDataFiles("./G1D1sgvsG1D1ngc")
G1D3sgvsG1D3ngc_G1D3wcvsG1D3ngc <- readDataFiles("./G1D3sgvsG1D3ngc")
G1D5sgvsG1D5ngc_G1D5wcvsG1D5ngc <- readDataFiles("./G1D5sgvsG1D5ngc")
G1D7sgvsG1D7ngc_G1D7wcvsG1D7ngc <- readDataFiles("./G1D7sgvsG1D7ngc")


UPonly <- list(DAG1 = G1D1sgvsG1D1ngc_G1D1wcvsG1D1ngc[["G1D1sgvsG1D1ngc_DEG_up.tsv"]]$gene_id,
               DAG3 = G1D3sgvsG1D3ngc_G1D3wcvsG1D3ngc[["G1D3sgvsG1D3ngc_DEG_up.tsv"]]$gene_id,
               DAG5 = G1D5sgvsG1D5ngc_G1D5wcvsG1D5ngc[["G1D5sgvsG1D5ngc_DEG_up.tsv"]]$gene_id,
               DAG7 = G1D7sgvsG1D7ngc_G1D7wcvsG1D7ngc[["G1D7sgvsG1D7ngc_DEG_up.tsv"]]$gene_id)


#### THIS GENERATES THE UPSET PLOT:

upset(fromList(UPonly),
      keep.order = T, 
      sets = c("DAG7", "DAG5", "DAG3", "DAG1"),
      mainbar.y.max = 700,
      sets.x.label = "No. differentially expressed genes",
      mainbar.y.label = "No. DEGs in common",
      #order = "freq",
      empty.intersections = "on",
      set_size.show = TRUE,
      set_size.scale_max = 1100,
      #number.angles = 90,
      text.scale = c(1.5, 1.5, 1.5, 1.5, 1.5, 1.5),
      point.size = 3.5, line.size = 2)



####  YOU CAN ALSO PLOT THE SAME DATA AS A VENN DIAGRAM:
library(VennDiagram)

# Helper function to display Venn diagram
display_venn <- function(UPonly, ...){
  
  grid.newpage()
  venn_object <- venn.diagram(UPonly, filename = NULL, ...)
  grid.draw(venn_object)
}

display_venn(
  UPonly,
  category.names = c("1 DAG" , "3 DAG" , "5 DAG", "7 DAG"),
  fill = c("#5BA8A0", "#CBE54E", "#94B447", "#5D6E1E")
)



######## Fig. 2b
# LETS EXTRACT grafted VS non-grafted DOWN-regulated genes only
######## 

G1D1sgvsG1D1ngc_G1D1wcvsG1D1ngc <- readDataFiles("./G1D1sgvsG1D1ngc")
G1D3sgvsG1D3ngc_G1D3wcvsG1D3ngc <- readDataFiles("./G1D3sgvsG1D3ngc")
G1D5sgvsG1D5ngc_G1D5wcvsG1D5ngc <- readDataFiles("./G1D5sgvsG1D5ngc")
G1D7sgvsG1D7ngc_G1D7wcvsG1D7ngc <- readDataFiles("./G1D7sgvsG1D7ngc")


DOWNonly <- list(DAG1 = G1D1sgvsG1D1ngc_G1D1wcvsG1D1ngc[["G1D1sgvsG1D1ngc_DEG_down.tsv"]]$gene_id,
                 DAG3 = G1D3sgvsG1D3ngc_G1D3wcvsG1D3ngc[["G1D3sgvsG1D3ngc_DEG_down.tsv"]]$gene_id,
                 DAG5 = G1D5sgvsG1D5ngc_G1D5wcvsG1D5ngc[["G1D5sgvsG1D5ngc_DEG_down.tsv"]]$gene_id,
                 DAG7 = G1D7sgvsG1D7ngc_G1D7wcvsG1D7ngc[["G1D7sgvsG1D7ngc_DEG_down.tsv"]]$gene_id)



upset(fromList(DOWNonly),
      keep.order = T, 
      sets = c("DAG7", "DAG5", "DAG3", "DAG1"),
      mainbar.y.max = 700,
      sets.x.label = "No. differentially expressed genes",
      mainbar.y.label = "No. DEGs in common",
      #order = "freq",
      empty.intersections = "on",
      set_size.show = TRUE,
      set_size.scale_max = 1100,
      #number.angles = 90,
      text.scale = c(1.5, 1.5, 1.5, 1.5, 1.5, 1.5),
      point.size = 3.5, line.size = 2)


####  YOU CAN ALSO PLOT THE SAME DATA AS A VENN DIAGRAM:
####  library(VennDiagram)

# Helper function to display Venn diagram
display_venn <- function(DOWNonly, ...){
  
  grid.newpage()
  venn_object <- venn.diagram(DOWNonly, filename = NULL, ...)
  grid.draw(venn_object)
}

display_venn(
  DOWNonly,
  category.names = c("1 DAG" , "3 DAG" , "5 DAG", "7 DAG"),
  fill = c("#5BA8A0", "#CBE54E", "#94B447", "#5D6E1E")
)







# ####
# #### THE FOLLOWING CODE MODIFIES THE UPSETR PACKAGE TO MOVE THE SET SIZE TO THE RIGHT OF THE PLOT
# #### UNCOMMENT with "shift/command/c ON A MAC IN RSTUDIO

# ######################
# library(grid)
# library(gridExtra)
# ######################
# 
# NoAttBasePlot <- function (legend, size_plot_height, Main_bar_plot, Matrix_plot, 
#                            hratios, Size_plot, query_legend, set_metadata, set_metadata_plots, 
#                            newpage) {
#   top <- 1
#   bottom <- 100
#   if ((!is.null(legend)) && (query_legend != tolower("none"))) {
#     if (query_legend == tolower("top")) {
#       top <- 3
#       bottom <- 102
#       legend_top <- 1
#       legend_bottom <- 3
#       size_plot_height <- (size_plot_height + 2)
#     }
#     else if (query_legend == tolower("bottom")) {
#       legend_top <- 101
#       legend_bottom <- 103
#     }
#   }
#   if (is.null(set_metadata)) {
#     matrix_and_mainbar_right <- 100
#     matrix_and_mainbar_left <- 21
#     size_bar_right <- 20
#     size_bar_left <- 1
#   }
#   else if (!is.null(set_metadata)) {
#     matrix_and_mainbar_right <- set_metadata$ncols + 100
#     matrix_and_mainbar_left <- set_metadata$ncols + 21
#     size_bar_right <- set_metadata$ncols + 20
#     size_bar_left <- set_metadata$ncols + 1
#     metadata_right <- set_metadata$ncols
#     metadata_left <- 1
#   }
#   if (newpage) {
#     grid::grid.newpage()
#   }
#   if ((!is.null(legend)) && (query_legend != tolower("none"))) {
#     if (query_legend == tolower("top")) {
#       pushViewport(viewport(layout = grid.layout(102, matrix_and_mainbar_right)))
#     }
#     else if (query_legend == tolower("bottom")) {
#       pushViewport(viewport(layout = grid.layout(103, matrix_and_mainbar_right)))
#     }
#   }
#   else if ((is.null(legend)) || (query_legend == tolower("none"))) {
#     pushViewport(viewport(layout = grid.layout(100, matrix_and_mainbar_right)))
#   }
#   # Modified
#   vp = UpSetR:::vplayout(top:bottom, 1:(matrix_and_mainbar_right-matrix_and_mainbar_left))
#   pushViewport(vp)
#   grid.draw(arrangeGrob(Main_bar_plot, Matrix_plot, heights = hratios))
#   popViewport()
#   # Modified
#   vp = UpSetR:::vplayout(size_plot_height:bottom, (matrix_and_mainbar_right-matrix_and_mainbar_left-1):96)
#   pushViewport(vp)
#   grid.draw(arrangeGrob(Size_plot))
#   popViewport()
#   if (!is.null(set_metadata)) {
#     for (i in 1:length(set_metadata_plots)) {
#       if (i != 1) {
#         metadata_left <- 1 + metadata_right
#         metadata_right <- metadata_right + set_metadata$plots[[i]]$assign
#       }
#       else {
#         metadata_left <- 1
#         metadata_right <- set_metadata$plots[[i]]$assign
#       }
#       vp = UpSetR:::vplayout(size_plot_height:bottom, metadata_left:metadata_right)
#       pushViewport(vp)
#       grid.draw(arrangeGrob(set_metadata_plots[[i]]))
#       popViewport()
#     }
#   }
#   if ((!is.null(legend)) && (query_legend != tolower("none"))) {
#     vp = UpSetR:::vplayout(legend_top:legend_bottom, matrix_and_mainbar_left:matrix_and_mainbar_right)
#     pushViewport(vp)
#     grid.draw(arrangeGrob(legend))
#     popViewport()
#   }
# }
# 
# Make_size_plot <- function (Set_size_data, sbar_color, ratios, ylabel, scale_sets, 
#                             text_scale, set_size_angle, set_size.show, set_size.scale_max, 
#                             set_size.number_size) {
#   if (length(text_scale) > 1 && length(text_scale) <= 6) {
#     x_axis_title_scale <- text_scale[3]
#     x_axis_tick_label_scale <- text_scale[4]
#   }
#   else {
#     x_axis_title_scale <- text_scale
#     x_axis_tick_label_scale <- text_scale
#   }
#   if (ylabel == "Set Size" && scale_sets != "identity") {
#     ylabel <- paste("Set Size", paste0("( ", 
#                                        scale_sets, " )"))
#     if (scale_sets == "log2") {
#       Set_size_data$y <- log2(Set_size_data$y)
#     }
#     if (scale_sets == "log10") {
#       Set_size_data$y <- log10(Set_size_data$y)
#     }
#   }
#   if (!is.null(set_size.number_size)) {
#     num.size <- (set_size.number_size/2.845276) * x_axis_tick_label_scale
#   }
#   else {
#     num.size <- (7/2.845276) * x_axis_tick_label_scale
#   }
#   Size_plot <- (ggplot(data = Set_size_data, aes_string(x = "x", 
#                                                         y = "y")) + geom_bar(stat = "identity", colour = sbar_color, 
#                                                                              width = 0.4, fill = sbar_color, position = "identity") + 
#                   scale_x_continuous(limits = c(0.5, (nrow(Set_size_data) + 
#                                                         0.5)), breaks = c(0, max(Set_size_data)), expand = c(0, 
#                                                                                                              0)) + theme(panel.background = element_rect(fill = "white"), 
#                                                                                                                          plot.margin = unit(c(-0.11, -1.3, 0.5, 0.5), "lines"), 
#                                                                                                                          axis.title.x = element_text(size = 8.3 * x_axis_title_scale), 
#                                                                                                                          axis.text.x = element_text(size = 7 * x_axis_tick_label_scale, 
#                                                                                                                                                     vjust = 1, hjust = 0.5), axis.line = element_line(colour = "gray0"), 
#                                                                                                                          axis.line.y = element_blank(), axis.line.x = element_line(colour = "gray0", 
#                                                                                                                                                                                    size = 0.3), axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
#                                                                                                                          panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + 
#                   xlab(NULL) + ylab(ylabel) + coord_flip())
#   if (set_size.show == TRUE) {
#     Size_plot <- (Size_plot + geom_text(aes(label = y, vjust = 0.5, 
#                                             hjust = 1.2, angle = set_size_angle), size = num.size))
#   }
#   if (scale_sets == "log10") {
#     if (!is.null(set_size.scale_max)) {
#       Size_plot <- (Size_plot + scale_y_continuous(limits = c(set_size.scale_max, 
#                                                               0), trans = log10_reverse_trans()))
#     }
#     else {
#       Size_plot <- (Size_plot + scale_y_continuous(trans = log10_reverse_trans()))
#     }
#   }
#   else if (scale_sets == "log2") {
#     if (!is.null(set_size.scale_max)) {
#       Size_plot <- (Size_plot + scale_y_continuous(limits = c(set_size.scale_max, 
#                                                               0), trans = log2_reverse_trans()))
#     }
#     else {
#       Size_plot <- (Size_plot + scale_y_continuous(trans = log2_reverse_trans()))
#     }
#   }
#   else {
#     if (!is.null(set_size.scale_max)) {
#       Size_plot <- (Size_plot + scale_y_continuous(limits = c(set_size.scale_max, 
#                                                               0), trans = "reverse"))
#     }
#     else {
#       # Modified
#       #Size_plot <- (Size_plot + scale_y_continuous(trans = "reverse"))
#     }
#   }
#   Size_plot <- ggplot_gtable(ggplot_build(Size_plot))
#   return(Size_plot)
# }
# 
# assignInNamespace(x="NoAttBasePlot", value=NoAttBasePlot, ns="UpSetR")
# assignInNamespace(x="Make_size_plot", value=Make_size_plot, ns="UpSetR")

#############################################




