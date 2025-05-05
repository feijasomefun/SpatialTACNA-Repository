###

# Recaculate CNB for all 273 samples 
# With TACNA profiles 383 TCs 
# Alterantive is in script 6

# Adjust TCs included in CNB 
# 1 TC high in normal removed , 9531 
# 2 TCs present in normal (bulk) tissue ; 9531, TC 562 is Y and therefore not removed(sig difference men and women)

# Small gene sets TCs were removed using sensitivity cut off for max (gene) value of 80 and a ratio (between first and second max) of 2.5  - 
# Genes with common names in more than 30% of max gene and in less than 50% of rest were also removed - 

# Create TACNA profiles with TCs < abs z score 3 are removed before creating profiles
# Save these profiles 
# Shift the TACNA profiles so the minimum score - smth becomes 0 
# Save this version 

# CNB is based on absolute colmeans 


###

# load packages 
library(data.table)
library(Seurat)
library(SeuratObject)
library(ggplot2)

# Directories 
INPUT_dir <- "/scratch/hb-bioinfo-fehrmann/Feija_Files/Stage_2/GSEA_finalize_CNA_TCs/Part_7_CNB_383_for_all_273_samples/Input_dir"
INPUT_CNB_dir <- "/scratch/hb-bioinfo-fehrmann/Feija_Files/cna_burdon/CNA_burden_RDS_output/Copy_number_burden_save/CNA_burden_RDS"

OUTPUT_dir <- "/scratch/hb-bioinfo-fehrmann/Feija_Files/Stage_2/GSEA_finalize_CNA_TCs/Part_7_CNB_383_for_all_273_samples/Output_dir"
DEGR_input_dir <- "/scratch/hb-bioinfo-fehrmann/Feija_Files/DEGR"

###

# 1 Load files 

###

# Import EC RDS files 
projected_samples_d <- file.path(INPUT_CNB_dir)
projected_samples_dir <- list.files(projected_samples_d, recursive = TRUE, full.names = TRUE)
print(projected_samples_dir)


# Update CNA TCs included
# Import CNA burdon ICA 
CNA_TC_ICA_dir <- file.path(INPUT_dir,"CNA_burdon_TC_ICA")
CNA_TC_ICA_f <- data.frame(fread(CNA_TC_ICA_dir), row.names = 1)

# change TC names 
colnames(CNA_TC_ICA_f) <- paste0("TC", sub("consensus.independent.component.", "", colnames(CNA_TC_ICA_f)))


# Copy number burden size factors doesnt work anymore 
# GSMS to plot and their pt.size.factor 
# GSM7211257 <- pt.size.factor = 2.3
# GSM7782699 <- pt.size.factor = 2
# GSM7850822 <- pt.size.factor = 2.1
# CC1 <- pt.size.factor = 1.7
# CRC2 <- pt.size.factor = 1.8
# EC8 <- pt.size.factor = 1.4 
# GSM5420752 <- pt.size.factor = 1.6
# GSM5494478 <- pt.size.factor = 2
# GSM5833536 <- pt.size.factor = 1.4
# GSM6506115 <- pt.size.factor = 1.9
# GSM7498812 <- pt.size.factor = 1.3
# GSM7990477 <- pt.size.factor = 1.4
# LC2 <- pt.size.factor = 1.7
# M1 <- pt.size.factor = 1.4

pt.size.map <- c(
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

###

# 2 Calculate CNB for each ST sample 

###


# Run the one without the TACNA 

f_generate_CNB <- function(
    projected_samples_dir,
    CNA_TC_ICA,
    CNA_TCs_included_names, 
    use_TACNA_profile_for_CNB, # True or false 
    OUTPUT_dir,
    output_name,
    save_as_RDS){

  for (sample_dir in projected_samples_dir){
    
    ### 
    
    # 1 Format Sample RDS 
    
    ###
    
    # insert sample directory 
    #sample_dir <- projected_samples_dir[8]
    print(sample_dir)
    
    # identify sample GSM or id 
    GSM_id <- sub("_.*", "", basename(sample_dir))
    
    # read RDS file 
    visium = readRDS(sample_dir)
    
    # extract the matrix with the TC activity scores
    sample_TC <- visium@assays[["Spatial"]]@layers[["counts"]]
    sample_TC_mtx <- as.matrix(sample_TC)
    
    # retrieve row (feature) and column (cell) names
    TC_names <- visium@assays[["Spatial"]]@features
    spot_names <- visium@assays[["Spatial"]]@cells
    length(TC_names) # should be 549 for all TCs 
    length(spot_names)
    dim(sample_TC_mtx)
    
    # gives spots (in columns) and TCs (should be 549 in rows) their names 
    colnames(sample_TC_mtx) <- rownames(spot_names)
    rownames(sample_TC_mtx) <- rownames(TC_names)
    
    # check sample 
    is.numeric(sample_TC_mtx)
    class(sample_TC_mtx)
    
    
    
    ###
    
    # 2 Re-calculate CNB 
    
    ###
    
    
    # Only use included TCs - automatically removes previous CNB 
    sample_TC_mtx_s <- sample_TC_mtx[which(rownames(sample_TC_mtx) %in% CNA_TCs_included_names),]
    
    print(length(rownames(sample_TC_mtx_s))) # should have 383
    
    # Set TCs with absolute activity scores less than 3 to 0 
    sample_TC_mtx_s_truncated = ifelse(abs(sample_TC_mtx_s)>3, sample_TC_mtx_s, 0)
    sample_TC_mtx_s_truncated = sample_TC_mtx_s
    
    # Ensure CNA_TC_ICA only has selected TCs and is a matrix
    CNA_TC_subset <- as.matrix(CNA_TC_ICA[, which(colnames(CNA_TC_ICA) %in% CNA_TCs_included_names)])
    
    if (use_TACNA_profile_for_CNB == TRUE){
      # matrix multiplication CNA burdon ICA and sample 
      sample_CNA_burdon_s <- CNA_TC_subset %*% sample_TC_mtx_s_truncated[colnames(CNA_TC_subset),]
    } else {sample_CNA_burdon_s <- sample_TC_mtx_s_truncated[colnames(CNA_TC_subset),]}
    
    # store the TACNA profiles
    # tacna_name <- paste0(GSM_id,"_TACNA_profile.txt")
    # write.table(sample_CNA_burdon_s, file = file.path(OUTPUT_dir,"EC_8_383_TCs_TACNA_profiles",tacna_name), sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
    
    
    # Absolute mean average of TACNA profiles 
    # use absolute mean value (for both positive and negative)
    average_CNA_burden_s <- colMeans(abs(sample_CNA_burdon_s))
    # create CNB scale so its from 0 (min) to 1 (max)
    CNB_383 = (average_CNA_burden_s - min(average_CNA_burden_s))/(max(average_CNA_burden_s) - min(average_CNA_burden_s))
    
    # Check the spots are still alligned 
    CNB_383 <- t(data.frame(CNB_383))
    identical(colnames(CNB_383),colnames(sample_TC_mtx))
    
    ###
    
    # 3 Create New RDS and plot 
    
    ###
    
    # Check the spatial spots have the same order 
    
    # New counts with burdon in last row # use old complete TC matrix 
    count_with_burdon_s <- rbind(sample_TC_mtx, CNB_383)
    
    # identical(count_with_burdon_c,count_with_burdon_s)
    
    # Automatize spatial spot size based on sample 
    
    # if (GSM_id == "EC1") { 
    #   spot_size <- 38
    # } else if (GSM_id == "EC2" || GSM_id == "EC3") { 
    #   spot_size <- 40
    # } else if (GSM_id == "EC5") { 
    #   spot_size <- 85
    # } else if (GSM_id == "EC6" || GSM_id == "EC7" || GSM_id == "EC8" || GSM_id == "EC9"|| GSM_id == "EC4") {
    #   spot_size <- 80  # for all other EC samples 
    # } else  {spot_size <- 1300}  # for all other EC samples
    # #}
    
    
    new.seurat.object = CreateSeuratObject(counts = count_with_burdon_s, assay = "Spatial" )
    new.seurat.object@images$image = new(
      
      Class = 'VisiumV1'
      
      ,assay = "spatial"
      
      ,key = "image_"
      
      ,coordinates = visium@images[["image"]]@coordinates
      
      ,image = visium@images[["image"]]@image
      
      ,scale.factors = visium@images[["image"]]@scale.factors
      
      ,spot.radius =  visium@images[["image"]]@spot.radius # spot_radius
      
    )
    
    # Change spot size 
    new.seurat.object@images$image@scale.factors$spot = visium@images$image@scale.factors$fiducial
    
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
    #custom_colors_p <- colorRampPalette(c("#005494", "#1a87d9", "#abd9e9", "#e4ffd1","#ffffbf", "#ffca96", "#f7421e", "#d7191c","#870002"))(length(breaks_profile))
    custom_colors_p <- colorRampPalette(c("#005494", "#1580d1", "#a2def2", "#e4ffd1","#fffac9", "#ffd7b0","#ff7559", "#d7191c","#870002"))(length(breaks_profile))
    
    
    pt.size.factor <- pt.size.map[[GSM_id]]
    if (is.null(pt.size.factor)) pt.size.factor <- 1.4 # default if not found
    
    # default_colors <- c("white", "blue", "red")  # Seurat uses a variant of magma
    CNA_Burden_plot <- SpatialFeaturePlot(new.seurat.object
                                          ,slot = "counts"
                                          ,features = "CNB-383"
                                          ,pt.size.factor = pt.size.factor) + #,pt.size.factor = spot_size)+
    scale_fill_gradientn(colors = custom_colors_p)  # Set limits manually
    
    
    sample_image <- SpatialFeaturePlot(new.seurat.object
                                       ,slot = "counts"
                                       ,features = "CNB-383"
                                       ,pt.size.factor = 0) 
    
    sample_image <- sample_image +NoLegend()
    
    # Plot 
    #sample_image
    #CNA_Burden_plot
    combined_plot <- sample_image + CNA_Burden_plot
    print(combined_plot)
    
    # Create directory 
    # Create the new output folder inside OUTPUT_dir
    new_output_dir <- file.path(OUTPUT_dir, output_name)
    if (!dir.exists(new_output_dir)) {
      dir.create(new_output_dir, recursive = TRUE)  # Ensure parent directories exist
    }
    
    # Create plot 
    basename <- paste0(GSM_id,"_",output_name)
    basename_pdf <- paste0(basename,".pdf")
    
    filename <- file.path(new_output_dir,basename_pdf)
    ggsave(
      filename = filename,  # file name for the saved image
      plot = combined_plot,            # the Seurat/ggplot object
      width = 14, height = 10
    )
    # 
    # # Export seurat object with feature plot 
    if (save_as_RDS == TRUE){
    saveRDS(new.seurat.object, file = file.path(new_output_dir, paste0(basename,".rds")))} #replace with rds ID but check spacing works
    
  }

}



# Subset for finalized TC 383

CNA_TCs_included <- data.frame(fread(file.path(INPUT_dir, "final_cna_TC_TC9531_removed_v2.tsv")),row.names = TRUE)
CNA_TCs_included_383_names <- rownames(CNA_TCs_included)
length(CNA_TCs_included_383_names) 

# Get updated ICA dataframe 
CNA_TC_ICA_383 <- CNA_TC_ICA_f[,which(colnames(CNA_TC_ICA_f)%in% CNA_TCs_included_383_names)]

f_generate_CNB(projected_samples_dir,
               CNA_TC_ICA = CNA_TC_ICA_383,
               CNA_TCs_included_names = CNA_TCs_included_383_names , 
               use_TACNA_profile_for_CNB = TRUE, # True or false
               OUTPUT_dir, 
               output_name = "TC_383_CNB_WITH_TACNA",
               save_as_RDS = TRUE) # true or false)

# Create selection of 14 high CNB samples 
selected_14_samples_folder_dir <- file.path(INPUT_dir,"selected_ST_samples_high_CNB_variance")
selected_14_samples_dir <- list.files(selected_14_samples_folder_dir, recursive = TRUE, full.names = TRUE)
selected_14_output_dir <- file.path(OUTPUT_dir,"CNB_14_groups")

f_generate_CNB(projected_samples_dir = selected_14_samples_dir,
               CNA_TC_ICA = CNA_TC_ICA_383,
               CNA_TCs_included_names = CNA_TCs_included_383_names ,
               use_TACNA_profile_for_CNB = TRUE, # True or false
               OUTPUT_dir = selected_14_output_dir,
               output_name = "TC_383_CNB_WITH_TACNA_14",
               save_as_RDS = TRUE) # true or false)


# Move all RDS files

# list rds files 
projected_rds_files <- list.files(OUTPUT_dir,full.names = TRUE, pattern = ".rds", recursive = TRUE)
dir_rds_folder <- file.path(OUTPUT_dir,"CNA_burden_RDS/")
#dir.create(dir_rds_folder)

# move each rds file to new folder
sapply(projected_rds_files, function(rds) {
  file.rename(rds, file.path(dir_rds_folder, basename(rds)))
})


# Move all TC plots 

# list pdf files 
projected_TC_plots <- list.files(OUTPUT_dir,full.names = TRUE, pattern = ".pdf",recursive = TRUE)
dir_TC_plots <- file.path(OUTPUT_dir,"CNA_burden_plots/")

#dir.create(dir_TC_plots)
# move each TC plot to new folder 
sapply(projected_TC_plots, function(plot) {
  file.rename(plot, file.path(dir_TC_plots, basename(plot)))
})


