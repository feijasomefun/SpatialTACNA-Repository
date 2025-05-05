###

# Create neat PDF plot for all ECs for any given feature 
# WARNING CNB 383 is called average-CNA-burden-s 
# CNB 549 is called average-CNA-burden 
# will redo and rename to 383 to double check 


###
# load packages 
library(data.table)
library(Seurat)
library(SeuratObject)
library(ggplot2)
library(patchwork)

# Directories 
INPUT_dir <- "/scratch/hb-bioinfo-fehrmann/Feija_Files/Stage_2/GSEA_finalize_CNA_TCs/part_6_TC_group_comparison/Input_dir"
OUTPUT_dir <- "/scratch/hb-bioinfo-fehrmann/Feija_Files/Stage_2/GSEA_finalize_CNA_TCs/part_6_TC_group_comparison/Output_dir"


###

# 1 Load files 

###

# Load in final 383 TCs 
final_383_TCs <- fread(file.path(INPUT_dir,"final_cna_TC_TC9531_removed_v2.tsv"))
names_383_TCs <- final_383_TCs$Component
length(names_383_TCs) # should be 383 

# Load 8 EC RDS files with 383 CNB 
# import folder with RDS files 
projected_samples_d <- file.path(INPUT_dir,"TC_383_update")
projected_samples_dir <- list.files(projected_samples_d, recursive = TRUE, full.names = TRUE, pattern = ".rds")
print(projected_samples_dir)


###

# Function to plot all samples in one neat pdf plot 

###


f_plot_samples_together <- function(projected_samples_dir,TC_to_plot,OUTPUT_dir,file_name){
  
  
  # Define breaks
  breaks_profile <- c(
    seq(-6, -5, length.out=2), # dark to mid
    seq(-5, -3, length.out=5),  # mid to light blue
    seq(-3, -1.5, length.out=4),  # light blue light green

    seq(-1.5, 0, length.out=3), # light green to yellow

    seq(0, 1.5, length.out=3), # yellow to orange

    seq(1.5, 3, length.out=4),  # orange to light red
    seq(3, 5, length.out=5),    # red to bright red
    seq(5, 6, length.out=2)) # bright to dark
  breaks_profile <- sort(unique(breaks_profile))  # Ensure uniqueness

  # Generate colors
  custom_colors_p <- colorRampPalette(c("#005494", "#1580d1", "#a2def2", "#e4ffd1","#fffac9", "#ffd7b0","#ff7559", "#d7191c","#870002"))(length(breaks_profile))

  # Create plot list
  plot_list <- list()
  
  for (sample_dir in projected_samples_dir){
    
    ### 
    
    # 1 Format Sample RDS 
    
    ###
    
    # insert sample directory 
    #sample_dir <- projected_samples_dir[8]
    print(sample_dir)
    
    # identify sample GSM or id 
    #GSM_id <- substr(basename(sample_dir),1,10)
    GSM_id <- sub("_.*", "", basename(sample_dir))
    
    # read RDS file 
    visium = readRDS(sample_dir)
    
    
    # Automatize spatial spot size based on sample 
    if (GSM_id == "EC1") { 
      spot_size <- 38
    } else if (GSM_id == "EC2" || GSM_id == "EC3") { 
      spot_size <- 40
    } else if (GSM_id == "EC5" || GSM_id == "EC7") { 
      spot_size <- 85
    } else {
      spot_size <- 80  # for all other EC samples 
    }
    
    
    ###
    
    # Create plot 
    
    if (TC_to_plot == "CNB-383"){
      
      TC_plot <- SpatialFeaturePlot(visium,
                                    slot = "counts",
                                    features = TC_to_plot,
                                    pt.size.factor = spot_size) +
        scale_fill_gradientn(colors = custom_colors_p,
                             oob = scales::squish,  # <-- This clips values to limits
                             values = rescale(breaks_profile, to = c(0, 1)),
                             limits = c(0, 1))  +
        theme(aspect.ratio = 1) + 
        coord_fixed() +
        ggtitle(GSM_id) +
        theme(plot.title = element_text(hjust = 0)) 
      
      
    } else if (TC_to_plot == "average-CNA-burden"){
      
      TC_plot <- SpatialFeaturePlot(visium,
                                    slot = "counts",
                                    features = TC_to_plot,
                                    pt.size.factor = spot_size) +
        scale_fill_gradientn(colors = custom_colors_p)  +
        theme(aspect.ratio = 1) + 
        coord_fixed() +
        ggtitle(GSM_id) +
        theme(plot.title = element_text(hjust = 0))} else{
    
    TC_plot <- SpatialFeaturePlot(visium,
                                    slot = "counts",
                                    features = TC_to_plot,
                                    pt.size.factor = spot_size) +
      scale_fill_gradientn(colors = custom_colors_p,
                           oob = scales::squish,  # <-- This clips values to limits
                           #values = rescale(breaks_profile, to = c(0, 1)),
                           limits = c(-6, 6),
                           labels = c("-6 <", "-3","0","3", "> 6"))  +
      theme(aspect.ratio = 1) + 
      coord_fixed() +
      ggtitle(GSM_id) +
      theme(plot.title = element_text(hjust = 0))}
    
    # Set limits manually 
    
    # Plot all ECs together 
    sample_image <- SpatialFeaturePlot(visium
                                       ,slot = "counts"
                                       ,features = "average-CNA-burden"
                                       ,pt.size.factor = 0) +NoLegend()
    
    TC_plot_with_image <- sample_image | TC_plot+
      plot_annotation(title = GSM_id) &
      theme(
        plot.title = element_text(hjust = -1, vjust = -10),
        plot.title.position = "plot"
      )
    print(TC_plot_with_image)
    
    plot_list[[GSM_id]] <- TC_plot_with_image
    
  }

  combined_plots <- wrap_plots(plotlist = plot_list, ncol = 2)
  
  # Create plot
  basename <- file_name
  
  filename <- file.path(OUTPUT_dir,basename)
  ggsave(
    filename = filename,  # file name for the saved image
    plot = combined_plots,            # the Seurat/ggplot object
    width = 10, height = 14
  )
  
  # Export seurat object with feature plot 
  # saveRDS(new.seurat.object, file = file.path(OUTPUT_dir, paste0(GSM_id,"c_TC","_fig2.rds"))) #replace with rds ID but check spacing works
}


###

# 2 Plot example TCs for fig 2 

# TCs to show 
# TC373 - good for contrasting groups
# TC2875 - good for CNB/tumour regions 
# TC1446 - small number of genes in TC 
# TC648 - good for lack of presence 

# TCs high in EC1 
# TC7844
# TC530 
# TC3794
# TC3591


f_plot_samples_together(
  projected_samples_dir = projected_samples_dir[1:8], 
  TC_to_plot =  c("TC373"),
  OUTPUT_dir = OUTPUT_dir, 
  file_name = "TC373_All_EC1_8_v2.pdf"
)


f_plot_samples_together(
  projected_samples_dir = projected_samples_dir[1:8], 
  TC_to_plot =  c("TC2875"),
  OUTPUT_dir = OUTPUT_dir, 
  file_name = "TC2875_All_EC1_8_v2.pdf"
)

# Note average CNB s is updated --> i should change the name to smth not confusing


f_plot_samples_together(
  projected_samples_dir = projected_samples_dir[1:8], 
  TC_to_plot =  c("CNB-383"),
  OUTPUT_dir = OUTPUT_dir, 
  file_name = "TC_383_CNB_All_EC1_8_v2.pdf"
)

# Note average CNB s is updated --> i should change the name to smth not confusing 
f_plot_samples_together(
  projected_samples_dir = projected_samples_dir[1:8], 
  TC_to_plot =  c("average-CNA-burden"),
  OUTPUT_dir = OUTPUT_dir, 
  file_name = "TC_549_CNB_All_EC1_8_v2.pdf"
)

