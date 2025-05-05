## Version 2 remove unnecessary addition of number of samples 


f_create_output_folder <- function(name_samples_folder)
{
  
  # Automatize the naming of the rds file (remove spaces)
  id_name <- gsub(" ", "_", name_samples_folder)
  id_name <- gsub(")", "", id_name)
  print(id_name)
  
  # create results folder with sample name 
  
  # check results directory
  # remove brackets from file name - avoid download problems 
  dir_results <- here("Results")
  print(dir_results)
  
  # create new folder in results directory
  dir.create(file.path(dir_results, id_name))
  
  qc_folder <- "qc_reports"
  # create new folder for qc reports 
  dir.create(file.path(dir_results, id_name,qc_folder))
  
  # create output directory 
  OUTPUT_DIR = here("Results",id_name)
  print(OUTPUT_DIR)
  
  return(OUTPUT_DIR)
}