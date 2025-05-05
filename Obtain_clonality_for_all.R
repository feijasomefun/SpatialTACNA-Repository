library(HiClimR)
library(Seurat)
library(data.table)
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("ConsensusClusterPlus")

library(ConsensusClusterPlus)
image_to_analyse = "EC5"


for(image_to_analyse in paste0("EC", 2:9))
{
  time2 = proc.time()[3]
  dir.create(file.path("/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Students/Feija Somefun internship/Manuscript/Main_figures/Fig_4_clonality/Results/", image_to_analyse), showWarnings = F)
  
  current_rds_file = readRDS(paste0('/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Students/Feija Somefun internship/Manuscript/Main_figures/Fig_2_spatial/Fig2_panel_C_CNB_with_tumour_region/Fig_EC_tumour_region/Input_data/',image_to_analyse,'_TC_383_CNB_WITH_TACNA.rds'))
  
  mix_mat = as.matrix(current_rds_file@assays$Spatial$counts)
  dim(mix_mat)
  
  degr_info = data.frame(fread('/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Students/Feija Somefun internship/Manuscript/Main_figures/Fig_2_spatial/Fig2_panel_C_CNB_with_tumour_region/Fig_EC_tumour_region/Input_data/DEGR_cut_383.tsv'), row.names = 1)
  length(unique(degr_info$Component))
  
  mix_mat = mix_mat[unique(degr_info$Component),]
  dim(mix_mat)
  
  mix_mat_cor = fastCor(mix_mat, nSplit = 25)
  
  
  setwd(paste0('/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Students/Feija Somefun internship/Manuscript/Main_figures/Fig_4_clonality/Results/', image_to_analyse, "/"))
  # Assuming final_cor_ordered is your data matrix
  
  time1 = proc.time()[3]
  results <- ConsensusClusterPlus(
    as.dist(1-mix_mat_cor),
    maxK = 10, # Maximum number of clusters to evaluate
    reps = 100, # Number of resampling iterations
    pItem = 0.8, # Proportion of items to sample in each iteration
    pFeature = 1, # Proportion of features to sample in each iteration
    title = paste("Consensus Clustering",image_to_analyse), # Title for output files
    innerLinkage = "ward.D2",
    finalLinkage = "ward.D2",
    clusterAlg = "hc", # Clustering algorithm (hierarchical clustering)
    distance = "pearson", # Distance metric
    seed = 1234, # Seed for reproducibility
    plot = "png", # Output plot format
    writeTable = FALSE
    # plot = NULL # Output plot format
  )
  
  print((proc.time()[3]-time1)/60)
  
  
  # Step 1: Compute area under CDF for each K
  get_cdf_area <- function(consensus_matrix) {
    # Extract only the lower triangle (excluding diagonal)
    vals <- consensus_matrix[lower.tri(consensus_matrix)]
    
    # Compute empirical CDF
    cdf_vals <- ecdf(vals)
    
    # Approximate area under the CDF curve (integration over [0,1])
    x <- seq(0, 1, length.out = 2000)
    y <- cdf_vals(x)
    auc <- trapz(x, y)  # Use trapezoidal integration
    
    return(auc)
  }
  
  # Trapezoidal integration function
  trapz <- function(x, y) {
    sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
  }
  
  # Step 2: Loop through your results and compute areas
  cdf_areas <- sapply(2:length(results), function(k) {
    get_cdf_area(results[[k]]$consensusMatrix)
  })
  
  # Step 3: Compute delta area
  delta_area <- diff(cdf_areas)
  delta_delta_area <- diff(delta_area)
  chosen_number_of_clusters = min(which(delta_delta_area>0))+3
  
  clones = results[[chosen_number_of_clusters]]$consensusClass 
  mix_mat = rbind(mix_mat, clones)
  
  
  new.seurat.object = CreateSeuratObject(counts = mix_mat, assay = "Spatial" )
  new.seurat.object@images$image = new(
    
    Class = 'VisiumV1'
    
    ,assay = "spatial"
    
    ,key = "image_"
    
    ,coordinates = current_rds_file@images[["image"]]@coordinates
    
    ,image = current_rds_file@images[["image"]]@image
    
    ,scale.factors = current_rds_file@images[["image"]]@scale.factors
    
    ,spot.radius =  current_rds_file@images[["image"]]@spot.radius # spot_radius
    
  )
  
  # Change spot size 
  new.seurat.object@images$image@scale.factors$spot = current_rds_file@images$image@scale.factors$fiducial
  
  
  clones = t(clones)
  rownames(clones) = "clones"
  write.table(clones,  file = paste0('/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Students/Feija Somefun internship/Manuscript/Main_figures/Fig_4_clonality/Results/', image_to_analyse, "/",image_to_analyse,"_clones_per_spot.txt"), sep = "\t",quote = FALSE)
  
  pdf(paste0('/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Students/Feija Somefun internship/Manuscript/Main_figures/Fig_4_clonality/Results/', image_to_analyse, "/",image_to_analyse,"_clones_per_spot.pdf"))
  print(SpatialFeaturePlot(new.seurat.object
                     ,slot = "counts"
                     ,features = "clones"))
  dev.off()
  
  print(image_to_analyse)
  print((proc.time()[3]-time2)/60)
  
}



