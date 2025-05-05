### explore 10 clusters in one sample 


# load packages 
library(data.table)
library(Seurat)
library(SeuratObject)
library(ggplot2)
library(patchwork)


### plot clonality 

INPUT_dir <- "/scratch/hb-bioinfo-fehrmann/Feija_Files/Stage_2/GSEA_finalize_CNA_TCs/Part_13_clonality"
Output_dir <- "/scratch/hb-bioinfo-fehrmann/Feija_Files/Stage_2/GSEA_finalize_CNA_TCs/Part_13_clonality/Output_data"



# Sample of interest GSM6592136

# Load 8 EC RDS files with 383 CNB 
# import folder with RDS files 
selected_samples <- file.path(INPUT_dir,"cluster_images_updated","GSM6592136_clusters_per_spot_2_to_10_clusters.rds") 


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
  
  

  # insert sample directory 
  sample_dir <- selected_samples[1]
  print(sample_dir)
  
  # identify sample GSM or id 
  #GSM_id <- substr(basename(sample_dir),1,10)
  GSM_id <- sub("_.*", "", basename(sample_dir))
  
  # read RDS file 
  visium = readRDS(sample_dir)
  
  # For first plot 
  #first_plot = TRUE
  plot_list <- list()
  
  for (current_feature in TC_to_plot){
    
    ### 
    
    # 1 Format Sample RDS 
    
    ###
    
    
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
    
    if (current_feature == "CNB-383"){
      
      TC_plot <- SpatialFeaturePlot(visium,
                                    slot = "counts",
                                    features = current_feature,
                                    pt.size.factor = pt.size.factor) +
        scale_fill_gradientn(colors = custom_colors_p,
                             oob = scales::squish,  # <-- This clips values to limits
                             values = rescale(breaks_profile, to = c(0, 1)),
                             limits = c(0, 1))  +
        theme(aspect.ratio = 1) + 
        coord_fixed() +
        ggtitle(GSM_id) +
        theme(plot.title = element_text(hjust = 0)) 
      
      
    } else{
      
      TC_plot <- SpatialFeaturePlot(visium,
                                    slot = "counts",
                                    features = current_feature,
                                    pt.size.factor = pt.size.factor) +
        #scale_fill_gradientn(colors = c("#078cb0","#fff27d","#ff8a5c"))  +
        theme(aspect.ratio = 1) + 
        coord_fixed() +
        ggtitle(GSM_id) +
        #NoLegend()+
        theme(plot.title = element_text(hjust = +2))} 
    TC_plot
    # Set limits manually 
    
    # Plot all ECs together 
    sample_image <- SpatialFeaturePlot(visium
                                       ,slot = "counts"
                                       ,features = "TC8"
                                       ,pt.size.factor = 0) +NoLegend()
    
    TC_plot_with_image <- TC_plot+
      plot_annotation(title = GSM_id) &
      theme(
        plot.title = element_text(hjust = -1, vjust = 0),
        plot.title.position = "plot"
      )
    print(TC_plot_with_image)
    
    plot_list[[current_feature]] <- TC_plot
    
  }
  
  combined_plots <- wrap_plots(plotlist = plot_list, ncol = 5)
  
  # Create plot
  basename <- file_name
  
  filename <- file.path(OUTPUT_dir,basename)
  ggsave(
    filename = filename,  # file name for the saved image
    plot = combined_plots,            # the Seurat/ggplot object
    width = 10, height = 4
  )
  
  # Export seurat object with feature plot 
  # saveRDS(new.seurat.object, file = file.path(OUTPUT_dir, paste0(GSM_id,"c_TC","_fig2.rds"))) #replace with rds ID but check spacing works
}



clusters_to_plot <- c("ClusterSet2clusters","ClusterSet3clusters",  "ClusterSet4clusters",  
                      "ClusterSet5clusters",  "ClusterSet6clusters",
                      "ClusterSet7clusters" , "ClusterSet8clusters" , "ClusterSet9clusters",  "ClusterSet10clusters")


# First half 
f_plot_samples_together(
  projected_samples_dir = selected_samples, 
  TC_to_plot =  clusters_to_plot,
  OUTPUT_dir = Output_dir, 
  file_name = "2_10_clusters.pdf"
)

