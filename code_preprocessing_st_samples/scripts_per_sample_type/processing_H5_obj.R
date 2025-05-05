library(Matrix)
library(Seurat)
library(data.table)
library(here)
library(R.utils)
library(knitr)
library(rmarkdown)
library(rhdf5)
source("https://raw.githubusercontent.com/loipf/coding_snippets/master/R/small_functions.R")

####

#1: INSERT samples FOLDER NAME and directory   
# FOLDERCEPTION based of updated processing 3 and pathfinder function 
# image type selected is high 

####
name_samples_folder <- "29_nasopharyngael_cancer"

# set directory to samples folder
samples_dir = here("Data",name_samples_folder)

print(samples_dir)
sample_files = list.files(path = samples_dir, pattern = "GSM")
print(sample_files)

# Find sample ID 
samples_id = unique(sapply(sample_files,function(x){
  
  # splits each file name x at underscores ("_")
  # return the first element of the resulting list.
  y = strsplit(x, "_")[[1]]
  
  # in this case only GSM that had h5 object 
  z = paste(y[1], sep = "_")
}))
print(samples_id)


# Unzip gz files 
zipped_files_dir <- list.files(path = samples_dir, pattern = "\\.gz$", full.names = TRUE)
for (zip_file in zipped_files_dir) {
  gunzip(zip_file, remove = TRUE)
}

sample_files <- list.files(path = samples_dir, pattern = "GSM.*\\.tar$", full.names = TRUE)

# Updated unzip function 
# can take tar of unzip since they have the same arguments 
f_zip_extraction_dir <- here("Code", "f_flexible_zip_extraction_v2.R")
print(f_zip_extraction_dir)
source(f_zip_extraction_dir)
f_zip_extraction(sample_files,samples_dir, unzip_function_name = untar)


# update sample files so its only the names and not full directory
sample_files <- list.files(path = samples_dir, pattern = "GSM")

####

# 2: CHECK FILE TYPE
# number of samples, presence H5 object and folder flatness

####

# insert function for checking sample file characteristics
f_check_duplicates_dir <- here("Code", "f_check_duplicates.R")
print(f_check_duplicates_dir)
source(f_check_duplicates_dir)

f_file_stats <- f_check_duplicates(sample_files)
multiples_samples_status <- f_file_stats$multiple_samples_status
folder_flat_status <- f_file_stats$folder_flat
h5_status <- f_file_stats$h5_status
number_of_samples <- f_file_stats$amount_unique_samples

# CHECK: the number of samples should be the same as number of samples id
# used in create output but dont seem to use it - check to delete 
correct_sample_number <- number_of_samples == length(samples_id)
print(correct_sample_number)

#SUMMARY FILE STATUS
print(paste0("File status: multiple samples is ", multiples_samples_status,
             ", flat folders ", folder_flat_status, " and H5 status is ", h5_status))

###

#3: PREPARE files for pre-processing 

###

#INSERT location of the NCBI file 

# set directory to NCBI folder
NCBI_dir <- here("NCBI")
# add file name to folder directory  
NCBI_file_name <- "/hgnc_gene_anno_merged.tsv"
NCBI_file_directory <- paste0(NCBI_dir,NCBI_file_name)

# fread NCBI
gene_anno = fread(NCBI_file_directory)
# Ensure NCBI gene ID's are characters 
gene_anno$NCBI_id = as.character(gene_anno$NCBI_id)

####

#4: CREATE OUTPUT folder and directory 
# automatized based on the sample folder input name 

####

f_output_folder_dir = here("Code", "f_create_output_folder_v2.R")
print(f_output_folder_dir)
# source the function from directory 
source(f_output_folder_dir)
# run function
OUTPUT_DIR <- f_create_output_folder(name_samples_folder)


####

#5: INSERT names required for seurat object 
# unzip and rename these files

####

#H5ls (rhdf5 package)



# Naming Conventions 
# dont need to unzip 
name_tissue_image <- "_tissue_hires_image.png" # insert here
name_scale_factors <- "_scalefactors_json.json"
name_tissue_position <- "_tissue_positions_list.csv"

name_h5_object <- "_filtered_feature_bc_matrix.h5"

### 
#6: PROCESS SAMPLES into tsv files with PC1 removed and RDS files 
# insert directory RDS function 

# check it works after downloading files 
RDS_function_dir = here("Code", "f_process_H5_samples_into_RDS_5.R")
source(RDS_function_dir)

# set the input of the function 
f_remove_PC1_form_RDS(
  samples_id, 
  samples_dir,
  name_h5_object,
  name_tissue_image,
  name_scale_factors,
  name_tissue_position,
  OUTPUT_DIR,
  gene_anno
)

