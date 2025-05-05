library(data.table)

gsea_tcga_100ev = data.frame(fread("/Users/arkajyotibhattacharya/Library/CloudStorage/OneDrive-UMCG/Data/TCGA_ICA_100ExplainedVariance_100runs/GSEA/enrichment_matrix_c2_cp_reactome_v2023_1.tsv"), row.names = 1)

library(Seurat)
sample_id = paste0("EC",1)
current_sample = readRDS(paste0('/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Spatial transcriptomics/Results/EC_samples_seurat_objects_with_TCGA_100EV_components/TCGA_components_on_',sample_id,'.rds'))
mix_mat = as.matrix(current_sample@assays$Spatial$counts)
dim(mix_mat)

clones = data.frame(fread(paste0('/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Students/Feija Somefun internship/Manuscript/Main_figures/Fig_6_clonality/Results/',sample_id,'/',sample_id,'_clones_per_spot.txt')))
rownames(clones) = clones$V1
clones$V1 = NULL

clones = t(clones)
rownames(clones) = sapply(rownames(clones), function(x){gsub("\\.", "-",x)})
sum(rownames(clones)==colnames(mix_mat))

cnb_file = readRDS('/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Students/Feija Somefun internship/Manuscript/Main_figures/Fig_4_contrasting_TC_activity/Input_data/RDS_383_CNB/EC1_TC_383_CNB_WITH_TACNA.rds')
cnb_file = as.matrix(cnb_file@assays$Spatial$counts)
sum(rownames(clones)==colnames(cnb_file))

cnb = cnb_file["CNB-383",]

write.table(rownames(gsea_tcga_100ev), file = "/Users/arkajyotibhattacharya/Downloads/reactome_genesets.txt", sep = "\t", quote = FALSE)

immune_genesets = data.frame(fread('/Users/arkajyotibhattacharya/Downloads/immune_genesets.txt'))
gsea_tcga_100ev_immune_genesets = gsea_tcga_100ev[immune_genesets$V2,]

gsea_tcga_100ev_immune_genesets_max = apply(gsea_tcga_100ev_immune_genesets,2,max)

summary(gsea_tcga_100ev_immune_genesets_max)

gsea_tcga_100ev_immune_genesets_max_indicator = ifelse(abs(gsea_tcga_100ev_immune_genesets)>5.56,1,0)
gsea_tcga_100ev_immune_genesets_max_indicator_colsums = colSums(gsea_tcga_100ev_immune_genesets_max_indicator)
table(gsea_tcga_100ev_immune_genesets_max_indicator_colsums)
names(gsea_tcga_100ev_immune_genesets_max_indicator_colsums) = sapply(names(gsea_tcga_100ev_immune_genesets_max_indicator_colsums),function(x){gsub("consensus.independent.component.", "TC",x)})
names(gsea_tcga_100ev_immune_genesets_max_indicator_colsums)[which(gsea_tcga_100ev_immune_genesets_max_indicator_colsums >0)] 
rownames(mix_mat) = sapply(rownames(mix_mat),function(x){gsub("consensus.independent.component.", "TC",x)})
mix_mat_immune_TCs = mix_mat[names(gsea_tcga_100ev_immune_genesets_max_indicator_colsums)[which(gsea_tcga_100ev_immune_genesets_max_indicator_colsums >0)] ,]

clone_numbers = sort(unique(clones[,1]))

cor_mat = matrix(NA, length(clone_numbers), nrow(mix_mat_immune_TCs))

for(i in 1:length(clone_numbers))
{
  
}
