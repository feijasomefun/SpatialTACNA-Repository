### 

# 21-08-2024 
# the untar and unzip function can be used interchangeably 
# since they take the same arguments 
# so they were combined and the unzip function needs to be given

###


f_zip_extraction <- function(sample_files,samples_dir,unzip_function_name){
  
  # for each sample in sample set 
  for (sample_dir in sample_files){
    
    #sample_dir <- "/scratch/hb-bioinfo-fehrmann/Feija_Files/Data/1_ovarian_cancer/GSM6506117_SP8_spatial.zip"
    
    # Split sample dir at GSM and extract next 7 letters for new folder 
    gsm_sample_id <- regmatches(sample_dir, regexpr("GSMT\\d{4}", sample_dir))
    sample_exit_dir <- file.path(samples_dir,gsm_sample_id)
    dir.create(sample_exit_dir)
    
    # Unzip files 
    # unzip_function_name(sample_dir,exdir = sample_exit_dir)
    
    # list all files in new subfolder
    list_sample_files_dir <- list.files(sample_exit_dir, recursive = TRUE, full.names = TRUE)
    
    # Add GSM to front of file names 
    for(file_dir in list_sample_files_dir){
      # Save basename of file 
      file_basename <- basename(file_dir)
      
      # Add the GSM ID at the front of the file basename 
      new_file_name <- paste0(gsm_sample_id,"_",file_basename)
      
      # Files will be put in the main directory of the data files
      new_file_dir <- file.path(samples_dir,new_file_name)
      
      # Rename file
      file.rename(file_dir, new_file_dir)
    }
    
    # delete empty GSM folder
    unlink(sample_exit_dir, recursive = TRUE)
  }
}
