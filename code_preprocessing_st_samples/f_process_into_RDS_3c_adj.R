# BASE FUNCTION FOR GENERATING RDS FILES AND PC1 REMOVED 
f_remove_PC1_form_RDS <- function(
    samples_id, 
    samples_dir,
    name_tissue_image,
    name_scale_factors,
    name_tissue_position,
    name_barcodes,
    name_features,
    name_matrix,
    OUTPUT_DIR,
    gene_anno
){
  
  ### IMPORT required functions 
  f_pathfinder_dir = here("Code", "f_pathfinder_v2.R")
  source(f_pathfinder_dir)
  quality_control_dir = here("Code", "f_quality_control_2.R")
  source(quality_control_dir)
  
  first_sample <- TRUE
  for(current_sample_id in samples_id) {
    
    # File selection 
    #current_sample_id <- "GSM7211257"
    # Update sample files, GSM will be used as ID  
    sample_files <- list.files(samples_dir)
    
    ###
    
    # 1 PATHFINDER and CHECK essential files 
    
    ###
    
    # TISSUE POSITION
    tissue_position.path <- f_pathfinder(current_sample_id,filename = name_tissue_position,sample_files, samples_dir,return_name = FALSE)
    # mapping file is for knowing where each barcode (row and column) is located 
    
    mapping_file <- read.csv(tissue_position.path, header = FALSE, row.names = 1)
    
    # MAPPING FILE: header presence 
    #check first whether first row contains word tissue
    # this indicates headers should be true 
    if (grepl("tissue", mapping_file[1,1]) == TRUE){
      # header present is true if values are characters
      header_present <- is.character(mapping_file[1,1])
      paste(mapping_file[1,])
    } else{header_present <- FALSE}
    
    # if the first column in first row is a character there is a header so headers should be on 
    # if not there is no header and a value shouldn't be accidentally removed 
    mapping_file <- read.csv(tissue_position.path, header = header_present, row.names = 1)
    colnames(mapping_file) <- c("tissue",   "row" ,     "col"  ,    "imagerow", "imagecol")
    
    # MATRIX 
    matrix.path <- f_pathfinder(current_sample_id,filename = name_matrix,sample_files, samples_dir,return_name = FALSE)
    # extract raw counts from matrix file 
    
    # Assign column names as barcodes and features as rownames 
    if (name_barcodes == name_features & name_features == name_matrix){
      
      raw_counts_csv <- read.csv(matrix.path, header = TRUE, row.names =1)
      raw_counts <- as.matrix(raw_counts_csv)
      mat_filtered <- as.matrix(raw_counts_csv)
      
      # barcodes are columns
      barcode.names <- colnames(raw_counts)
      # features/genes are rows 
      feature.names <- rownames(raw_counts)
      # replace any dots with - in columnnames 
      colnames(raw_counts) <- gsub("\\.", "-", colnames(raw_counts))
      
    } else {
      
      raw_counts = as.matrix(readMM(file = matrix.path))
      mat_filtered =  as.matrix(readMM(file = matrix.path))
      
      barcode.path <- f_pathfinder(current_sample_id,filename = name_barcodes,sample_files, samples_dir,return_name = FALSE)
      features.path <- f_pathfinder(current_sample_id,filename = name_features,sample_files, samples_dir,return_name = FALSE)
      
      barcode.names <- read.delim(barcode.path, header = FALSE, stringsAsFactors = FALSE)
      feature.names <- read.delim(features.path, header = FALSE, stringsAsFactors = FALSE)
      
      # column names of raw count is barcodes' column V1 which are probe spot locations: ATTGTCGCAATACCTT-1
      colnames(raw_counts) <- barcode.names$V1
      # row names of raw count becomes features' column V2 which is (hugo) gene symbols
      rownames(raw_counts) <- feature.names$V2
    }
    
    # MATRIX checks 
    # check matrix dimensions 
    check_matrix_dim <- dim(mat_filtered)
    # check the sum of NAs in the matrix 
    check_matrix_na<- sum(is.na(mat_filtered))
    
    # RAW COUNTS checks 
    # check first 3 column names 
    check_column_names <- colnames(raw_counts)[1:3]
    # check the number of characters for column names (should be barcode with 18 characters)  
    nchar_colnames <- sapply(colnames(raw_counts), nchar)
    check_column_sum <- summary(nchar_colnames)
    raw_counts = raw_counts[names(which(table(rownames(raw_counts))==1)),]
    # check first 3 row names 
    check_row_names <- rownames(raw_counts)[1:3]
    # check number of characters in rownames(should be gene names)
    nchar_rownames <- sapply(rownames(raw_counts), nchar)
    check_row_sum <- summary(nchar_rownames)
    
    # check common samples: mapping file and tissue position
    # probe barcodes that are both in barcodes file (can change to matrix columns) and mapping file
    common_samples = intersect(rownames(mapping_file), colnames(raw_counts)) # change to colnames matrix 
    mapping_file = mapping_file[common_samples,]
    raw_counts_colsums = colSums(raw_counts)
    
    # select raw counts which that have a sum of more than 0 in columns
    # columns/barcodes that have all the same value are removed
    raw_counts = raw_counts[,names(which(raw_counts_colsums>0))]
    raw_counts_colsd = apply(raw_counts,2,sd)
    raw_counts = raw_counts[,names(which(raw_counts_colsd>0))]
    
    # Incase of removed columns ensure rows still match with mapping file
    # find common barcodes: rownames of mapping are barcodes colnames of rawcounts are barcodes 
    common_samples = intersect(rownames(mapping_file), colnames(raw_counts))
    mapping_file = mapping_file[common_samples,] # use mapping for RDS 
    raw_counts = as.matrix(raw_counts[,common_samples]) # genes (rows) x spot location (columns)
    
    # RAW COUNTS checks 
    check_dimensions_raw_count <- as.character(dim(raw_counts))
    # Extra checks:
    check_na_raw_count <- sum(is.na(raw_counts))
    check_row_dup_raw_count <- sum(duplicated(colnames(raw_counts)))
    
    ###
    
    #2 PERFORM PCA, remove PC1 export final ncbi matrix
    
    ###
    
    # Prepare raw counts for PCA : 
    # raw counts row sums should not be 0, these are removed 
    raw_counts_pca = raw_counts[rowSums(raw_counts)!=0,]
    # log normalize the raw counts: Feature counts for each cell are divided 
    # by the total counts for that cell and multiplied by the scale.factor. 
    # followed by natural log transformation
    #? what is the scale factor here
    # make sure a few very highly active genes are not overwhelming other signals 
    raw_counts_pca = Seurat::LogNormalize(raw_counts_pca)
    
    # Perform PCA
    pca_obj_princomp = princomp(raw_counts_pca, cor=T, fix_sign = F)  ### input: genes x samples
    # pca_obj_princomp$loadings == at6 eigenvector
    
    ### remove PC_1 with PCA 
    #pc1_loading_factor <- as.data.frame(pca_obj_princomp$loadings[,1])
    # check the sign of the PC1 loading factors
    # Loading factors is the activity of each PC at that spot
    freq_table_PC1 <- table(sign(pca_obj_princomp$loadings[,1])) # current sign positive or negative 
    
    # remove PC1: background effect 
    # matrix multiplication loading factors and transposed scores 
    # extract loadings except PC1, extact scores except PC1
    # multiply loading with transposed scores 
    # rows becomes spots and columns genes 
    raw_counts_pc1removed = pca_obj_princomp$loadings[,-1] %*% t(pca_obj_princomp$scores[,-1]) 
    
    # exports as tsv file 
    fwrite(as.data.frame(raw_counts_pc1removed), file.path(OUTPUT_DIR, paste0(current_sample_id,"_lognorm_pc1removed.tsv")), row.names=T, sep="\t")
    
    ### replace with entrez ncbi id
    # prevent overlapping/improper naming
    # match the colnames of raw counts PC1 removed with HGNC symbol
    # for those that match find equivelant NCBI ID 
    colnames_ncbi <- gene_anno$NCBI_id[match(colnames(raw_counts_pc1removed), gene_anno$HGNC_symbol)]
    # use NCBI ID as colnames 
    colnames(raw_counts_pc1removed) <- colnames_ncbi
    raw_counts_pc1removed <- raw_counts_pc1removed[, !is.na(colnames(raw_counts_pc1removed))]
    raw_counts_pc1removed <- raw_counts_pc1removed[, !duplicated(colnames(raw_counts_pc1removed))]
    print(paste0("final matrix ncbi: ", nrow(raw_counts_pc1removed)," samples x ",ncol(raw_counts_pc1removed)," genes"))
    
    # create raw counts file with PC1 removed 
    fwrite(as.data.frame(raw_counts_pc1removed), file.path(OUTPUT_DIR, paste0(current_sample_id,"_lognorm_pc1removed_ncbi.tsv")), row.names=T, sep="\t")
    
    ###
    
    # CHECK PC1 removed raw counts 
    # rows currently spots , genes currently cols 
    
    ###
    
    # SD for rows and columns should not be 0; show number that are 0 
    
    # rows are spots cols are genes 
    # perform SD of rows (1) 
    check_pc1_removed_sd_row <- apply(raw_counts_pc1removed,1,sd)
    # check length of named vector containing row SD of a 0 value 
    check_pc1_removed_sd_row <- length(which(check_pc1_removed_sd_row == 0))
    # perform SD of columns (2)
    check_pc1_removed_sd_col <- apply(raw_counts_pc1removed,2,sd)
    # check length of named vector containing col SD of a 0 value
    check_pc1_removed_sd_col <- length(which(check_pc1_removed_sd_col==0))
    
    # check  length of col with Sd of 0 
    check_pc1_removed_na <- sum(is.na(raw_counts_pc1removed))
    
    # Original names using pathfinder
    full_image.path <- f_pathfinder(current_sample_id,filename = name_tissue_image,sample_files, samples_dir, return_name = FALSE)
    full_image_name <- f_pathfinder(current_sample_id,filename = name_tissue_image,sample_files, samples_dir, return_name = TRUE)
    print(full_image_name)
    full_json_name <- f_pathfinder(current_sample_id,filename = name_scale_factors,sample_files, samples_dir, return_name = TRUE)
    print(full_json_name)
    full_json.path <- f_pathfinder(current_sample_id,filename = name_scale_factors,sample_files, samples_dir, return_name = FALSE)
    full_tissue_position_name <- f_pathfinder(current_sample_id,filename = name_tissue_position,sample_files, samples_dir, return_name = TRUE)
    print(full_tissue_position_name)
    
    # Check the correlation matrix of principle components 
    pc1_removed_scores <-as.data.frame(pca_obj_princomp$scores[,-1])
    # correlation matrix: activity of PCs should not be correlated 
    
    if (first_sample == TRUE) 
    {
      # number of PCs should be number of barcodes (-1 from PC1 removed)
      cormat = cor(pc1_removed_scores) #the middle ones should be 1 the rest 0
      # # correlation of 1 PC2-PC2 and the rest should be completely uncorrelated very near 0 
      cormat_sum <- sum(cormat)
      # # cormat sum, sum of number of barcodes 
      check_cormat_sum <- c(cormat_sum , nrow(raw_counts_pc1removed))
      first_sample <- FALSE
    } else{
      check_cormat_sum <- print("cormat was already run on the first sample")
    } 
    
    
    ###
    
    #3 Create, plot and export Seurat object
    
    ###
    
    # setwd so image selection works 
    setwd(samples_dir)
    
    # make temporary names for read10x function 
    file.rename(from = full_json_name, to = "scalefactors_json.json")
    # I convert any image to be called lowres 
    file.rename(from = full_image_name, to = "tissue_lowres_image.png")
    file.rename(from = full_tissue_position_name, to = "tissue_positions_list.csv")
    
    # spatial data goes wrong 
    # the mapping is incorrect 
    # create seurat object
    spatial_data = Read10X_Image(
      image.dir = file.path(samples_dir), 
      filter.matrix = FALSE, )
    
    # swap back to original name which has sample id
    file.rename(from = "tissue_lowres_image.png", to = full_image_name)
    file.rename(from = "scalefactors_json.json", to = full_json_name)
    file.rename(from = "tissue_positions_list.csv", to = full_tissue_position_name )
    
    
    # for hires original images 
    if(grepl("hires", full_image_name) == TRUE )
    {
      image <- png::readPNG(source = full_image.path)
      scale.factors <- jsonlite::fromJSON(txt = full_json.path)
      # change spot radius to be the radius of the hires scale factor 
      # this is done by reversing what is done in 10xgenomicsimages
      # times by low scale to remove low scale
      # next times by hi scale to get the hi scale unormalized radius 
      unnormalized.radius <- scale.factors$fiducial_diameter_fullres*scale.factors$tissue_hires_scalef
      spot.radius <- unnormalized.radius/max(dim(x = image))
      spatial_data@spot.radius = spot.radius
      # set lowres to be the hires scale factor to ensure correct image scaling
      spatial_data@scale.factors[["lowres"]] = spatial_data@scale.factors$hires
      
    }
    
    new.seurat.object = CreateSeuratObject(counts = raw_counts, assay = "Spatial" )
    new.seurat.object@images$image = new(
      Class = 'VisiumV1'
      ,assay = "spatial"
      ,key = "image_"
      ,coordinates = mapping_file
      ,image = spatial_data@image
      ,scale.factors = spatial_data@scale.factors
      ,spot.radius = spatial_data@spot.radius
    )
    
    processed_data <- list(
      # matrix filtered is raw counts without PC1 removed  
      mat_filtered = mat_filtered,
      feature_names = feature.names,
      barcode_names = barcode.names,
      spatial_data = new.seurat.object
    )
    
    # Save Spatial seurat object and files as RDS file for each sampleS 
    saveRDS(processed_data, file = file.path(OUTPUT_DIR, paste0(current_sample_id,".rds"))) #replace with rds ID but check spacing works
    
    # Plot housekeeping genes: GADPH, ACTB , B2M 
    gene_of_interest <- c("GADPH","ACTB","B2M")
    
    # Plot spatial feature - gene of interest to check image and scale factor 
    check_image_and_scale_plot <- SpatialFeaturePlot(new.seurat.object
                                                     ,slot = "counts"
                                                     , features = gene_of_interest)

    # Import quality control function 
    quality_control_dir = here("Code", "f_quality_control_2.R")
    source(quality_control_dir)
    
    
    quality_control <- f_quality_control_2(check_matrix_na, 
                                           # Raw Counts 
                                           check_column_names,nchar_colnames, check_column_sum,
                                           check_row_sum, check_row_names, 
                                           check_dimensions_raw_count, 
                                           # Raw counts check 
                                           check_na_raw_count,
                                           check_row_dup_raw_count,
                                           # PC1 checks 
                                           freq_table_PC1, check_pc1_removed_na, 
                                           check_pc1_removed_sd_row, check_pc1_removed_sd_col,
                                           # Check Correlation matrix PCs 
                                           check_cormat_sum)
    
    
    # Create QC directory
    qc_output_dir = paste0(OUTPUT_DIR,"/qc_reports")
    qc_report_dir <- here("Code","qc_report_template_2.Rmd")
    
    # Export QC to pdf file 
    generate_pdf_qc_report <- function(qc_output_dir, qc_report_dir, current_sample_id, quality_control,check_image_and_scale_plot) {
      output_file <- paste0(qc_output_dir,"/QC_report_", current_sample_id, ".pdf")
      rmarkdown::render(
        # create generic input rmarkdown file 
        input = qc_report_dir,
        output_file = output_file,
        # alter parameters of specific file 
        params = list(sample_name = current_sample_id, 
                      qc_result = quality_control,
                      qc_plot = check_image_and_scale_plot
        ),
        # create new environment from parent
        envir = new.env(parent = globalenv())
      )
    }
    
    generate_pdf_qc_report(qc_output_dir, qc_report_dir, current_sample_id, quality_control, check_image_and_scale_plot)
    
  }
}


# F quality control function - too hectic to include in processing 