# UPDATE: 11-06-2024 the folder depth wasnt working changed function to depend on duplications 
# flat folders have duplicates deep do not 
# insert function for checking sample file characteristics


f_check_duplicates <- function(sample_files){
  sample_list <- list()
  
  # check if files have .h5
  h5_status <- any(unlist(lapply(sample_files, function(file) grepl("\\.h5$", file))))
  #print(paste0("h5 status ", h5_status))
  
  # create list of files with only sample id 
  for (current_file in sample_files) {
    current_sample <- substr(current_file, 1, 10)
    sample_list <- append(current_sample, sample_list)
  }
  
  # check for duplicates 
  duplicates <- duplicated(sample_list)
  
  # check for unique duplicates of the same sample
  unique_duplicates <- unique(sample_list[duplicates])
  print(paste0("unique duplicates ",unique_duplicates))
  
  # size/length of list is number of unique samples 
  number_unique_duplicate_samples <-length(unique_duplicates)
  print(paste0("number unique duplicate samples ",number_unique_duplicate_samples))
  print(paste0("There are ", number_unique_duplicate_samples, " samples with duplications out of ", length(sample_list), " files"))
  
  # if there are no duplicates the folders are deep 
  if (number_unique_duplicate_samples == 0){
    folder_flat <- FALSE
  }else{
      folder_flat <- TRUE
  }
  
  # check number of unique GSM files 
  unique_files <- unique(sample_list)
  
  #count number of unique file names
  amount_unique_samples <- length(unique_files)
  print(paste0("number of unique files ", amount_unique_samples))
  #print(paste0("Multiple sample status is ", multiple_samples_status))
  multiple_samples_status <- sum(amount_unique_samples) > 1 
  
  # count list of duplicates
  num_duplicates <- sum(duplicates)
  print(paste0("There are ", num_duplicates," duplications in total"))
  
  
  # for every unique sample in unique duplicates 
  for (unique_sample in unique_duplicates){
    # count number of unique samples in sample list 
    count <- sum(sample_list == unique_sample)
    # if more than one folder of same sample ID the count is > 1
    # folderception 
    print(paste0("There are ", count, " files of the following sample: ", unique_sample)) 
  }
  
  output_list <- list(multiple_samples_status = multiple_samples_status, folder_flat = folder_flat, h5_status = h5_status, amount_unique_samples = amount_unique_samples)
  return(output_list)
}
