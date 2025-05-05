# updated pathfinder 19-08-2024 to seperate files 
# name_tissue_image, name_scale_factors, name_tissue_position have unzipped and zipped version
# need unzipped version for function but other files are still zipped 

f_pathfinder <- function(current_sample_id,filename,sample_files, samples_dir,return_name){
  
  #current_sample_id <- "GSM5494475"
  #filename <- name_tissue_position
  pattern <- paste0(current_sample_id, ".*", filename)
  
  # grepl file with this pattern
  sample_file_name <- sample_files[grepl(pattern, sample_files)]
  #print(sample_file_name)
  # since not all files have unzipped version 
  if (filename %in% c(name_tissue_image, name_scale_factors, name_tissue_position)){
    filtered_name <- sample_file_name[!grepl("\\.gz$", sample_file_name)]
  } else{
    filtered_name <- sample_file_name
  }
  #print(filtered_name)
  path <- file.path(samples_dir, filtered_name)
  if (return_name == TRUE){
    return(filtered_name)
  } else{
    return(path)
  }
}
