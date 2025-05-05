
###
# load packages 
library(data.table)
library(Seurat)
library(SeuratObject)
library(ggplot2)
library(patchwork)


### plot clonality 

INPUT_dir <- "/scratch/hb-bioinfo-fehrmann/Feija_Files/Stage_2/GSEA_finalize_CNA_TCs/Part_13_clonality"
Output_dir <- "/scratch/hb-bioinfo-fehrmann/Feija_Files/Stage_2/GSEA_finalize_CNA_TCs/Part_13_clonality/Output_data"


###

# 1 Load files 

###

# Load in final 383 TCs 

# Load 8 EC RDS files with 383 CNB 
# import folder with RDS files 
projected_samples_d <- file.path(INPUT_dir,"cluster_images_updated")
projected_samples_dir <- list.files(projected_samples_d, recursive = TRUE, full.names = TRUE, pattern = ".rds")
print(projected_samples_dir)
selected_samples <- projected_samples_dir[!grepl("^EC", basename(projected_samples_dir))]
print(selected_samples)



pt.size.map <- c(
  GSM6339638 = 1.6,
  GSM5420751 = 1.8,
  GSM6177612 = 1.8,
  GSM6177623 = 1.8,
  GSM7211257 = 2.3,
  GSM7782699 = 2,
  GSM7850822 = 2.1,
  CC1 = 1.7,
  CRC2 = 1.8,
  EC8 = 1.4,
  GSM5420752 = 1.6,
  GSM5494478 = 1.4,
  GSM5833536 = 1.4,
  GSM6506115 = 1.9,
  GSM7498812 = 1.3,
  GSM7990477 = 1.4,
  LC2 = 1.7,
  M1 = 1.4
)


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
  
  # Colors 
  custom_colors_p <- colorRampPalette(c("#005494", "#1580d1", "#a2def2", "#e4ffd1","#fffac9", "#ffd7b0","#ff7559", "#d7191c","#870002"))(length(breaks_profile))
  
  
  
  # For first plot 
  #first_plot = TRUE
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
    
    # Change spot size 
    # new.seurat.object@images$image@scale.factors$spot = visium@images$image@scale.factors$fiducial
    # Size factor still needs to be adjusted despite changed spot size
    if (GSM_id %in% names(pt.size.map)) {
      pt.size.factor <- pt.size.map[[GSM_id]]
    } else {
      pt.size.factor <- 1.4
    }
    ###
    
    # Create plot 
    
    if (TC_to_plot == "CNB-383"){
      
      TC_plot <- SpatialFeaturePlot(visium,
                                    slot = "counts",
                                    features = TC_to_plot,
                                    pt.size.factor = pt.size.factor) +
        scale_fill_gradientn(colors = custom_colors_p,
                             oob = scales::squish,  # <-- This clips values to limits
                             values = rescale(breaks_profile, to = c(0, 1)),
                             limits = c(0, 1))  +
        theme(aspect.ratio = 1) + 
        coord_fixed() +
        ggtitle(GSM_id) +
        theme(plot.title = element_text(hjust = 0)) 
      
      
    } else if (TC_to_plot == "ClusterSet3clusters"|| TC_to_plot == "clusters"){
      
      TC_plot <- SpatialFeaturePlot(visium,
                                    slot = "counts",
                                    features = TC_to_plot,
                                    pt.size.factor = pt.size.factor) +
        scale_fill_gradientn(colors = c("#078cb0","#fff27d","#ff8a5c"))  +
        theme(aspect.ratio = 1) + 
        coord_fixed() +
        ggtitle(GSM_id) +
        NoLegend()+
        theme(plot.title = element_text(hjust = +2))} else{
          
          TC_plot <- SpatialFeaturePlot(visium,
                                        slot = "counts",
                                        features = TC_to_plot,
                                        pt.size.factor = pt.size.factor) +
            scale_fill_gradientn(colors = custom_colors_p,
                                 oob = scales::squish,  # <-- This clips values to limits
                                 #values = rescale(breaks_profile, to = c(0, 1)),
                                 limits = c(-6, 6),
                                 labels = c("-6 <", "-3","0","3", "> 6"))  +
            theme(aspect.ratio = 1) + 
            coord_fixed() +
            ggtitle(GSM_id) +
            theme(plot.title = element_text(hjust = 0))}
    TC_plot
    # Set limits manually 
    
    # Plot all ECs together 
    sample_image <- SpatialFeaturePlot(visium
                                       ,slot = "counts"
                                       ,features = "TC8"
                                       ,pt.size.factor = 0) +NoLegend()
    
    TC_plot_with_image <- sample_image | TC_plot+
      plot_annotation(title = GSM_id) &
      theme(
        plot.title = element_text(hjust = -1, vjust = 0),
        plot.title.position = "plot"
      )
    print(TC_plot_with_image)
    
    plot_list[[GSM_id]] <- TC_plot_with_image
    
  }
  
  combined_plots <- wrap_plots(plotlist = plot_list, ncol = 2)
  
  # Create plot
  basename <- file_name
  height_to_plot <- length(projected_samples_dir)*2
  filename <- file.path(OUTPUT_dir,basename)
  ggsave(
    filename = filename,  # file name for the saved image
    plot = combined_plots,            # the Seurat/ggplot object
    width = 10, height = height_to_plot
  )
  
  # Export seurat object with feature plot 
  # saveRDS(new.seurat.object, file = file.path(OUTPUT_dir, paste0(GSM_id,"c_TC","_fig2.rds"))) #replace with rds ID but check spacing works
}


###

# 2 Plot example TCs for fig 2E


# First half 
f_plot_samples_together(
  projected_samples_dir = selected_samples[1:14], 
  TC_to_plot =  "ClusterSet3clusters",
  OUTPUT_dir = Output_dir, 
  file_name = "3_clusters_1_14.pdf"
)


# EC using clusterset 3 

selected_samples_EC <- projected_samples_dir[grepl("^EC", basename(projected_samples_dir))]

# First half 
f_plot_samples_together(
  projected_samples_dir = selected_samples_EC, 
  TC_to_plot =  "ClusterSet3clusters",
  OUTPUT_dir = Output_dir, 
  file_name = "3_clusters_EC1_8.pdf"
)


### EC samples plotting

projected_samples_EC_dir <- file.path(INPUT_dir,"cluster_images_updated")
projected_samples_EC <- list.files(projected_samples_EC_dir, recursive = TRUE, full.names = TRUE, pattern = ".rds")
print(projected_samples_EC)
# Make sure to get the cluster per spot one 
selected_samples <- projected_samples_EC[grepl("EC", basename(projected_samples_dir))&grepl("clusters_per_spot.rds", basename(projected_samples_dir))]
selected_samples_EC <- selected_samples[1:8]
print(selected_samples_EC)

# Second half 
f_plot_samples_together(
  projected_samples_dir = selected_samples_EC, 
  TC_to_plot =  "clusters",
  OUTPUT_dir = Output_dir, 
  file_name = "3_clusters_EC_1_8_perspot_file.pdf"
)



###

# Plot all 273 as a catalogue 

# Remove EC9 
projected_samples_dir_273 <- 	list.files(file.path(INPUT_dir, "cluster_images_all"), pattern = ".rds",full.names = TRUE)
projected_samples_dir_273 <- projected_samples_dir_273[!grepl("^EC", basename(projected_samples_dir_273))]




f_plot_all <- function(samples_per_page, feature_to_plot, projected_samples_dir, OUTPUT_dir){
  
  first_sample <- 1 
  last_sample <- first_sample + samples_per_page - 1
  
  
  while(first_sample <= length(projected_samples_dir)){
    
    file_name <- paste0(feature_to_plot, "_", first_sample, "_", last_sample,".pdf")
    
    f_plot_samples_together(
      projected_samples_dir = projected_samples_dir[first_sample:last_sample], 
      TC_to_plot = feature_to_plot,
      OUTPUT_dir = Output_dir, 
      file_name = file_name
    )
    
    # Update for next page
    first_sample <- last_sample + 1
    last_sample <- first_sample + samples_per_page - 1
    
    # To avoid out-of-bounds error:
    if (first_sample > length(projected_samples_dir)){
      break
    }
    if (last_sample > length(projected_samples_dir)){
      last_sample <- length(projected_samples_dir)
    }
  }
}

f_plot_all(samples_per_page = 8, 
           feature_to_plot = "clusters",
           projected_samples_dir = projected_samples_dir_273, 
           OUTPUT_dir = Output_dir)



