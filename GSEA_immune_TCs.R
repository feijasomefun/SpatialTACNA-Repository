library(data.table)
library(ggplot2)
library(cowplot)

gsea_tcga_100ev = data.frame(fread("/Users/arkajyotibhattacharya/Library/CloudStorage/OneDrive-UMCG/Data/TCGA_ICA_100ExplainedVariance_100runs/GSEA/enrichment_matrix_c2_cp_reactome_v2023_1.tsv"), row.names = 1)


immune_genesets = data.frame(fread('/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Students/Feija Somefun internship/Manuscript/Main_figures/Fig_6_clonality/Input data/immune_genesets.txt'))
gsea_tcga_100ev_immune_genesets = gsea_tcga_100ev[immune_genesets$V2,]

gsea_tcga_100ev_immune_genesets_max = apply(gsea_tcga_100ev_immune_genesets,2,max)

summary(gsea_tcga_100ev_immune_genesets_max)

gsea_tcga_100ev_immune_genesets_max_indicator = ifelse(abs(gsea_tcga_100ev_immune_genesets)>5.56,1,0)
gsea_tcga_100ev_immune_genesets_max_indicator_colsums = colSums(gsea_tcga_100ev_immune_genesets_max_indicator)
table(gsea_tcga_100ev_immune_genesets_max_indicator_colsums)
names(gsea_tcga_100ev_immune_genesets_max_indicator_colsums) = sapply(names(gsea_tcga_100ev_immune_genesets_max_indicator_colsums),function(x){gsub("consensus.independent.component.", "TC",x)})
colnames(gsea_tcga_100ev_immune_genesets) = sapply(colnames(gsea_tcga_100ev_immune_genesets),function(x){gsub("consensus.independent.component.", "TC",x)})

immune_TCs = names(gsea_tcga_100ev_immune_genesets_max_indicator_colsums)[which(gsea_tcga_100ev_immune_genesets_max_indicator_colsums >0)] 
gsea_tcga_100ev_immune_genesets_immunne_TCs = gsea_tcga_100ev_immune_genesets[,immune_TCs]

rownames(gsea_tcga_100ev_immune_genesets_immunne_TCs) = sapply(rownames(gsea_tcga_100ev_immune_genesets_immunne_TCs), function(x){gsub("REACTOME_", "", x)})
rownames(gsea_tcga_100ev_immune_genesets_immunne_TCs) = sapply(rownames(gsea_tcga_100ev_immune_genesets_immunne_TCs), function(x){strsplit(x , " -- ")[[1]][1]})
rownames(gsea_tcga_100ev_immune_genesets_immunne_TCs) = sapply(rownames(gsea_tcga_100ev_immune_genesets_immunne_TCs), function(x){gsub("_", " ", x)})
rownames(gsea_tcga_100ev_immune_genesets_immunne_TCs) = sapply(rownames(gsea_tcga_100ev_immune_genesets_immunne_TCs), function(x){substr(x, 1, 40)})

gsea_tcga_100ev_immune_genesets_immunne_TCs_rowSums= rowSums(gsea_tcga_100ev_immune_genesets_immunne_TCs)
gsea_tcga_100ev_immune_genesets_immunne_TCs = gsea_tcga_100ev_immune_genesets_immunne_TCs[which(gsea_tcga_100ev_immune_genesets_immunne_TCs_rowSums>0),]
pdf('/Users/arkajyotibhattacharya/Downloads/immune_TC_GSEA.pdf', width = 10, height = 12)
pheatmap::pheatmap(gsea_tcga_100ev_immune_genesets_immunne_TCs
                   ,clustering_distance_rows = "correlation"
                   ,clustering_distance_cols = "correlation"
                   ,clustering_method = "ward.D2"
)
dev.off()


