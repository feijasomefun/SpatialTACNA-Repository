
library(tidyverse)
library(stringr)
library(data.table)

###
# Input dir 

# We removed small gene sets TCs using sensitivity cut off for max value of 80 and a ratio,
# We removed genes with common names in more than 30% of max gene and in less than 50% of rest 
# Additionally 9531 still needs to be removed - its highly active in normal samples 

# Arko checked if genes overlap with z score genes in geset enrichment analyses
# only 9531 did which will be removed because it was high in normal 


###
INPUT_dir <- "/scratch/hb-bioinfo-fehrmann/Feija_Files/Stage_2/GSEA_finalize_CNA_TCs/Part_2_recalculate_CNB_EC_samples/Input_dir"
OUTPUT_dir <- "/scratch/hb-bioinfo-fehrmann/Feija_Files/Stage_2/GSEA_finalize_CNA_TCs/Part_5_finalizing_CNA_TCs/Exploration/Output_dir"
DEGR_input_dir <- "/scratch/hb-bioinfo-fehrmann/Feija_Files/DEGR"



# Load in CNA-TCs 
CNA_TCs_highest_z_gene_set <- data.frame(fread(file.path(INPUT_dir,"549_CNA_TCs_highest_score_reactome_GOBP.txt")), row.names = TRUE)
all_549_CNA_TCs <- rownames(CNA_TCs_highest_z_gene_set)


###

# 2 Experiment with first and second max 

###


# Load in DEGFR/TACNA results
# Anotation file which shows which TC belongs to which chromosome 
TC_position <- fread(file.path(DEGR_input_dir,"Genelevel_DEGR_TCGA_100_EV_All_Components_Amplified_or_Deleted_details_FDR_0.05_CL_0.5_state_deciding_cutoff_0.85_.txt"))
TC_position$Component <-  paste0("TC",TC_position$Component)
TC_position_cut <- TC_position[TC_position$Component %in% all_549_CNA_TCs,]


###

# 3 Remove common name / protein complex TCs 

###

# Step 1: extract common string 

extract_common_substring <- function(str) {
  # Use regex to capture the common substring; this will match words and numbers
  common_substring <- str_extract(str, "^([\\w\\s]+(?:[\\w]+))")  # Match words and numbers (no special chars)
  return(common_substring)
}


# Step 2: Apply the function to extract common substrings
TC_common_genes <- TC_position_cut %>%
  mutate(Common_Protein_Complex = sapply(GENENAME, extract_common_substring))


# Create new dataframe showing the max value 

TC_max_value_df <-  TC_position_cut %>% 
  group_by(Component) %>%  # Group by Component
  summarise(
    
    Max_Value = max(Value, na.rm = TRUE),
    Second_Max = max(Value[Value < max(Value, na.rm = TRUE)], na.rm = TRUE)
  )%>%
  mutate(Ratio = Max_Value / Second_Max)  # Calculate the ratio


# Step 3: function to find the longest common substring
longest_common_substring <- function(string1, string2) {
  n1 <- nchar(string1)
  n2 <- nchar(string2)
  
  # Create a matrix to store lengths of longest common suffixes of substrings
  matrix <- matrix(0, nrow = n1 + 1, ncol = n2 + 1)
  
  # Store the length and position of the longest common substring
  max_length <- 0
  lcs_end <- 0
  
  # Build the matrix
  for (i in 1:n1) {
    for (j in 1:n2) {
      if (substr(string1, i, i) == substr(string2, j, j)) {
        matrix[i + 1, j + 1] <- matrix[i, j] + 1
        if (matrix[i + 1, j + 1] > max_length) {
          max_length <- matrix[i + 1, j + 1]
          lcs_end <- i
        }
      }
    }
  }
  
  # If there is no common substring
  if (max_length == 0) {
    return(NULL)  # No common substring found
  }
  
  # Extract the longest common substring
  return(substr(string1, lcs_end - max_length + 1, lcs_end))
}

# Find the longest common substring
# lcs <- longest_common_substring(string1, string2)
# print(lcs)


unique_TCs = unique(TC_position_cut$Component)
TC_max_value_df$common_char_string_frequency = NA
TC_max_value_df$common_char_string = NA
TC_max_value_df$total_number_of_genes = NA
TC_max_value_df$number_of_characters_in_highest_gene_without_common_part = NA
TC_max_value_df$total_number_of_characters_in_highest_gene = NA

for(i in 1:length(unique_TCs))
{
  current_TC_detail = TC_position_cut[which(TC_position_cut$Component==unique_TCs[i]),]
  
  genenames = toupper(current_TC_detail$GENENAME)
  
  # common_substring_matrix = matrix(NA, length(genenames), length(genenames))
  # 
  # for(gene1 in 1:length(genenames))
  # {
  #   for(gene2 in gene1:length(genenames))
  #   {
  #     the_common_char = longest_common_substring(genenames[gene1], genenames[gene2])
  #     if(!is.null(the_common_char))
  #     {
  #       common_substring_matrix[gene1, gene2] = the_common_char
  #     }
  #     
  #   }  
  # }
  genenames_strsplit = sapply(genenames, function(x){strsplit(x, " ")[[1]]})
  highest_weighted_genename = genenames_strsplit[[which.max(current_TC_detail$Value)]]
  
  
  common_substring_array = array(NA, length(genenames))
  
  for(gene2 in c(1:length(genenames)))
  {
    the_common_char = intersect(highest_weighted_genename, genenames_strsplit[[gene2]])
    if(length(the_common_char)>0)
    {
      common_substring_array[gene2] = paste(the_common_char, collapse = " ")
    }
    
  }  
  
  TC_index = which(TC_max_value_df$Component==unique_TCs[i])
  common_substring_array = common_substring_array[which(common_substring_array!="")]
  TC_max_value_df$common_char_string_frequency[TC_index] = max(table(common_substring_array))+1
  TC_max_value_df$common_char_string[TC_index] = names(which.max(table(common_substring_array)))
  TC_max_value_df$total_number_of_genes[TC_index] = nrow(TC_position_cut[which(TC_position_cut$Component==unique_TCs[i])])
  TC_max_value_df$total_number_of_characters_in_highest_gene[TC_index] = nchar(paste(highest_weighted_genename, collapse = " "))
  TC_max_value_df$number_of_characters_in_highest_gene_without_common_part[TC_index] = TC_max_value_df$total_number_of_characters_in_highest_gene[TC_index] - nchar( TC_max_value_df$common_char_string[TC_index])
  print(i)
}


cut_off_ratios = c(11:100)/10
cut_off_max = c(2:10)*10

number_of_TCs_after_chopping_off = matrix(NA, 9,90)

for(i in 1:9)
{
  for(j in 1:81)
  {
    number_of_TCs_after_chopping_off[i,j] = length(intersect(which(TC_max_value_df$Max_Value > cut_off_max[i])
                                                             , which(TC_max_value_df$Ratio >cut_off_ratios[j] )))
    
  }
}

rownames(number_of_TCs_after_chopping_off) = paste("cut_off_max", cut_off_max)
colnames(number_of_TCs_after_chopping_off) = paste("cut_off_ratios", cut_off_ratios)


small_geneset_driven_TCs = TC_max_value_df$Component[intersect(which(TC_max_value_df$Max_Value >= 80)
                                                               , which(TC_max_value_df$Ratio >=2.5 ))]



rest_of_the_TCs =  TC_max_value_df$Component[which(! TC_max_value_df$Component%in%small_geneset_driven_TCs)]


# Remove TCs with a high number of characters 
TC_max_value_df_rest_of_the_TCs = TC_max_value_df[which(! TC_max_value_df$Component%in%small_geneset_driven_TCs),]
TC_max_value_df_rest_of_the_TCs$ratio_of_protein_complex_genes_in_evr = TC_max_value_df_rest_of_the_TCs$common_char_string_frequency/TC_max_value_df_rest_of_the_TCs$total_number_of_genes
TC_max_value_df_rest_of_the_TCs$ratio_of_non_common_characters = TC_max_value_df_rest_of_the_TCs$number_of_characters_in_highest_gene_without_common_part/TC_max_value_df_rest_of_the_TCs$total_number_of_characters_in_highest_gene

TC_max_value_df_rest_of_the_TCs_protein_complex = TC_max_value_df_rest_of_the_TCs[intersect(which(TC_max_value_df_rest_of_the_TCs$ratio_of_protein_complex_genes_in_evr>0.3)
                                                                                            , which(TC_max_value_df_rest_of_the_TCs$ratio_of_non_common_characters < 0.5)),]
                                                                                  


# Export table of protein complex 
write.table(TC_max_value_df_rest_of_the_TCs_protein_complex,file = file.path(OUTPUT_dir,"TC_max_value_df_rest_of_the_TCs_protein_complex.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)                                                                                

# There are 11 protein complexe TCs that are removed 
length(TC_max_value_df_rest_of_the_TCs_protein_complex)

# Export table of small TC driven - should be 154 
TC_max_value_small_geneset <- TC_max_value_df[TC_max_value_df$Component %in% small_geneset_driven_TCs,]

write.table(TC_max_value_small_geneset, file.path(OUTPUT_dir,"TC_max_value_small_geneset.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)


# Export finalized TC list 
final_cna_TC = TC_max_value_df_rest_of_the_TCs[which(!TC_max_value_df_rest_of_the_TCs$Component%in%TC_max_value_df_rest_of_the_TCs_protein_complex$Component),]
# Remove TC5931 
final_cna_TC_TC9531_removed <- final_cna_TC[which(final_cna_TC$Component != "TC9531"), ]

write.table(final_cna_TC_TC9531_removed, file.path(OUTPUT_dir, "final_cna_TC_TC9531_removed.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

###

# Plots to export 

###


plot(TC_max_value_df$total_number_of_genes, TC_max_value_df$Ratio, pch = 19, col = rgb(0,0,0,0.1))

# Sensitivity Plot # big jump 160 and 154 in the table 
plot(number_of_TCs_after_chopping_off[7, ])
# Export the sensitity table 
write.table(number_of_TCs_after_chopping_off, file.path(OUTPUT_dir, "number_of_TCs_after_chopping_off.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)


# Plot first and second max - so ratio 
second_Max_distribution <- ggplot(TC_max_value_df, aes(x = Max_Value, y = Second_Max)) +
  geom_point(color = "blue", size = 1, alpha = 0.4) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +  # Add reference line at Ratio = 1
  labs(title = "Scatter Plot of Ratios",
       x = "First max",
       y = "Second Max") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate labels for readability
second_Max_distribution



