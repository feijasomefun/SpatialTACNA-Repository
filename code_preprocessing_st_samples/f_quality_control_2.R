
f_quality_control_2 <- function(check_matrix_na, 
                               # Raw Count check 
                               check_column_names,nchar_colnames, check_column_sum,
                               check_row_sum, check_row_names, 
                               check_dimensions_raw_count, 
                               # Raw count check 
                               check_na_raw_count,
                               check_row_dup_raw_count,
                               # PC1 removed check
                               freq_table_PC1,check_pc1_removed_na, 
                               check_pc1_removed_sd_row, check_pc1_removed_sd_col,
                               # check matrix corelation of PCs 
                               check_cormat_sum
                               ){
  
  # Checks NAs in matrix file should be 0  
  p_check_matrix_na <- paste("The sum of NAs in matrix are:", check_matrix_na)
  
  ###
  
  # Check rawcounts table
  
  ###
  
  # Check columns 
  # Check the number of characters for column names which should be probe barcodes
  
  # Check names of first three columns in rawcounts
  p_check_RC_colnames <- check_column_names
  #p_check_RC_colnames <- c(paste("The names of the first three columns are:", paste(check_column_names, collapse = ", ")))
  
  # Automatize the checking of the length of column names - should be consistently 18 characters for barcodes 
  if (all(nchar_colnames == 18)){
    p_check_RC_len_colnames <- print("All the column names have 18 characters, which is the correct amount for barcodes")
  } else {
    p_check_RC_len_colnames <- print(" Error not all the column names are 18 which is abnormal for barcode length implying incorrect column names")
  }
  p_check_RC_colsum_raw <- summary(nchar_colnames)
  
  p_check_RC_dim <- c(paste("The dimension (rows * columns) of the raw counts table are: ", 
                                          paste(check_dimensions_raw_count,collapse = " * ")))
  
  # Rownames check
  # The number of characters in rownames should be gene names which have a higher max than the number of rows
  # Gene names can be very large much larger than for example 5 characteres in number of rows 30000 
  # This check indicates gene names are present instead of numerical 
  
  p_check_RC_row_sum_expl_1 <- paste(check_row_sum["Max."], "characters")
  p_check_RC_row_sum_expl_2 <- paste(check_dimensions_raw_count[1], "rows")
  p_check_RC_row_sum_expl <- paste(p_check_RC_row_sum_expl_1, p_check_RC_row_sum_expl_2, sep = " , ")
  #p_check_RC_row_sum_expl <- paste(check_row_sum["Max."], check_dimensions_raw_count[1], sep = " , ")
  
  p_check_RC_row_sum <- print(check_row_sum)
  
  # Check first names 
  p_check_RC_row_names <- print(paste("The names of the first three rows are:", paste(check_row_names, collapse = ", ")))
  
  
  # Checks length of frequency table: should only have one sign +1 or -1 
  if (length(freq_table_PC1) == 1) {
    sign_present <- names(freq_table_PC1)
    p_PC1_sign <- paste("The loading factor for PC1 are all of one sign, as shown by frequency count of +1 or -1:")
  } else {
    p_PC1_sign <- paste("The loading factor for PC1 have multiple signs, as shown by frequency counts of -1, +1:")
  }
  
  p_check_freq_PC1_sign <- c(p_PC1_sign ,freq_table_PC1)
  
  p_check_pc1_removed_na <- paste("The number of NA's in PC1 removed is:",check_pc1_removed_na)
  p_check_pc1_removed_sd_row <- paste("The number of rows with SD of 0 in PC1 removed is:",check_pc1_removed_sd_row)
  p_check_pc1_removed_sd_col <- paste("The number of columns with SD of 0 in PC1 removed is:",check_pc1_removed_sd_col)
  
  p_check_RC_na <- paste("The number of NAs in raw counts is:",check_na_raw_count)
  p_check_RC_row_dup <- paste("The number of duplicates in raw counts rows is:",check_row_dup_raw_count)
  
  
  # Correlation matrix sum 
  
  p_check_cormat_sum <- check_cormat_sum
  
  # Have output for checks 
  p_check_list <- list(p_check_matrix_na, 
                       # RAW counts (RC): 
                       p_check_RC_colnames,p_check_RC_len_colnames,p_check_RC_colsum_raw,
                       p_check_RC_row_sum_expl, p_check_RC_row_sum,p_check_RC_row_names,
                       # changed order 
                       p_check_RC_dim,p_check_RC_na,p_check_RC_row_dup,
                       # PC1 removed check 
                       p_check_freq_PC1_sign,p_check_pc1_removed_na, 
                       p_check_pc1_removed_sd_row, p_check_pc1_removed_sd_col, 
                       # Correlation matrix sum 
                       p_check_cormat_sum
                      )

  
  return(p_check_list)
}

