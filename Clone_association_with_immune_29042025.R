library(data.table)
library(ggplot2)
library(cowplot)
gsea_tcga_100ev = data.frame(fread("/Users/arkajyotibhattacharya/Library/CloudStorage/OneDrive-UMCG/Data/TCGA_ICA_100ExplainedVariance_100runs/GSEA/enrichment_matrix_c2_cp_reactome_v2023_1.tsv"), row.names = 1)

library(Seurat)
sample_id = paste0("EC",1)
current_sample = readRDS(paste0('/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Spatial transcriptomics/Results/EC_samples_seurat_objects_with_TCGA_100EV_components/TCGA_components_on_',sample_id,'.rds'))
mix_mat = as.matrix(current_sample@assays$Spatial$counts)
dim(mix_mat)
# '/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Students/Feija Somefun internship/Manuscript/Main_figures/Fig_6_clonality/Input data/cluster_images/EC2_clusters_per_spot.txt'
clones = data.frame(fread(paste0('/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Students/Feija Somefun internship/Manuscript/Main_figures/Fig_6_clonality/Input data/cluster_images/',sample_id,'_clusters_per_spot.txt')))
rownames(clones) = clones$V1
clones$V1 = NULL

clones = t(clones)
rownames(clones) = sapply(rownames(clones), function(x){gsub("\\.", "-",x)})
sum(rownames(clones)==colnames(mix_mat))

cnb_file = readRDS(paste0('/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Students/Feija Somefun internship/Manuscript/Main_figures/Fig_4_contrasting_TC_activity/Input_data/RDS_383_CNB/',sample_id,'_TC_383_CNB_WITH_TACNA.rds'))
cnb_file = as.matrix(cnb_file@assays$Spatial$counts)
sum(rownames(clones)==colnames(cnb_file))

cnb = cnb_file["CNB-383",]

# write.table(rownames(gsea_tcga_100ev), file = "/Users/arkajyotibhattacharya/Downloads/reactome_genesets.txt", sep = "\t", quote = FALSE)

immune_genesets = data.frame(fread('/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Students/Feija Somefun internship/Manuscript/Main_figures/Fig_6_clonality/Input data/immune_genesets.txt'))
gsea_tcga_100ev_immune_genesets = gsea_tcga_100ev[immune_genesets$V2,]

gsea_tcga_100ev_immune_genesets_max = apply(gsea_tcga_100ev_immune_genesets,2,max)

summary(gsea_tcga_100ev_immune_genesets_max)

gsea_tcga_100ev_immune_genesets_max_indicator = ifelse(abs(gsea_tcga_100ev_immune_genesets)>5.56,1,0)
gsea_tcga_100ev_immune_genesets_max_indicator_colsums = colSums(gsea_tcga_100ev_immune_genesets_max_indicator)
table(gsea_tcga_100ev_immune_genesets_max_indicator_colsums)
names(gsea_tcga_100ev_immune_genesets_max_indicator_colsums) = sapply(names(gsea_tcga_100ev_immune_genesets_max_indicator_colsums),function(x){gsub("consensus.independent.component.", "TC",x)})
immune_TCs = names(gsea_tcga_100ev_immune_genesets_max_indicator_colsums)[which(gsea_tcga_100ev_immune_genesets_max_indicator_colsums >0)] 
rownames(mix_mat) = sapply(rownames(mix_mat),function(x){gsub("consensus.independent.component.", "TC",x)})
mix_mat_immune_TCs = mix_mat[names(gsea_tcga_100ev_immune_genesets_max_indicator_colsums)[which(gsea_tcga_100ev_immune_genesets_max_indicator_colsums >0)] ,]
mix_mat_immune_TCs_with_cnb = as.data.frame(rbind(mix_mat_immune_TCs, cnb))
clone_numbers = sort(unique(clones[,1]))

cor_mat = matrix(NA, length(clone_numbers), nrow(mix_mat_immune_TCs))
colnames(cor_mat) = immune_TCs
rownames(cor_mat) = paste0("cluster",c(1:3))

cor_mat_spots = mix_mat_immune_TCs

for(i in 1:length(clone_numbers))
{
  spot_in_current_cluster = rownames(clones)[which(clones[,1]==i)]
  mix_mat_immune_TCs_current = mix_mat_immune_TCs_with_cnb[,spot_in_current_cluster]
  
  cor_mat_current= cor(t(mix_mat_immune_TCs_current))
  
  for(j in 1:dim(cor_mat_spots)[1])
  {
    cor_mat_spots[j,spot_in_current_cluster] = cor_mat_current["cnb", rownames(cor_mat_spots)[j]]
  }
  cor_mat[i,] = cor_mat_current["cnb",-length(cor_mat_current["cnb",])]
  
}
rownames(cor_mat_spots) = paste0(rownames(cor_mat_spots), "CorWithCnb")

mix_mat_with_cnb = rbind(mix_mat, cnb)
mix_mat_with_cnb_cor_with_immune_TCs = rbind(mix_mat_with_cnb, cor_mat_spots)
current_rds_file = readRDS(paste0('/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Students/Feija Somefun internship/Manuscript/Main_figures/Fig_4_contrasting_TC_activity/Input_data/RDS_383_CNB/',sample_id,'_TC_383_CNB_WITH_TACNA.rds'))

new.seurat.object = CreateSeuratObject(counts = mix_mat_with_cnb_cor_with_immune_TCs, assay = "Spatial")
new.seurat.object@images$image = new(
  Class = 'VisiumV1',
  assay = "spatial",
  key = "image_",
  coordinates = current_rds_file@images[["image"]]@coordinates,
  image = current_rds_file@images[["image"]]@image,
  scale.factors = current_rds_file@images[["image"]]@scale.factors,
  spot.radius = current_rds_file@images[["image"]]@spot.radius
)

new.seurat.object@images$image@scale.factors$spot = current_rds_file@images$image@scale.factors$fiducial
saveRDS(new.seurat.object, paste0('/scratch/hb-bioinfo-fehrmann/Feija_Files/Results/Clonality_all/cluster_images_all', "/", image_to_analyse, "_clusters_per_spot.rds"))
SpatialFeaturePlot(new.seurat.object, slot = "counts", features = "TC8708CorWithCnb", pt.size.factor = 1.6)

images_to_plot = sample_id
tcs_to_plot = c('raw',"cnb", immune_TCs,'raw',"cnb",rownames(cor_mat_spots))
plot_list_tcs = lapply(images_to_plot, function(curr_img){
  
  print(paste0("### start ", curr_img))
  
  plot_list_tcs_curr_tc = lapply(tcs_to_plot, function(curr_tc){
    
    curr_tc_nchar = nchar(curr_tc)
    
    if(curr_tc =='raw') {
      
      p = SpatialFeaturePlot(new.seurat.object, slot = "counts", features = "cnb", pt.size.factor = 1.6, alpha = 0)
      p = p + theme(plot.margin=margin(l=-0.5,unit="in"))
      
    } else if(curr_tc_nchar< 10){
      p = SpatialFeaturePlot(new.seurat.object, slot = "counts", features = curr_tc, pt.size.factor = 1.6)
      p = p + theme(plot.margin=margin(l=-0.5,unit="in"))
      
    }else {
      p = SpatialFeaturePlot(new.seurat.object, slot = "counts", features = curr_tc, pt.size.factor = 1.6)
      p = p + theme(plot.margin=margin(l=-0.5,unit="in"))
    }
    print(paste0("### ", curr_tc))
    
    
    p
    
  })
  
  plot_list_tcs_curr_tc
  
})

plot_list_tcs = unlist(plot_list_tcs, recursive = F)



pdf(paste0('/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Students/Feija Somefun internship/Manuscript/Main_figures/Fig_6_clonality/Plots/', sample_id,'_cnb_immune_association.pdf'), width = 3*length(tcs_to_plot)*2, height = 2*length(images_to_plot)*20)

plot_grid(plotlist = plot_list_tcs
          # , nrow=length(images_to_plot)
          , nrow=2
          , rel_widths = c(0.6,rep(1,length(tcs_to_plot)+2)), scale = 1)

dev.off()

pdf(paste0('/scratch/hb-bioinfo-fehrmann/Feija_Files/Results/Clonality_all/cluster_images_all', "/", image_to_analyse, "_clusters_per_spot.pdf"))
p1 = SpatialFeaturePlot(new.seurat.object, slot = "counts", features = "clusters", pt.size.factor = pt_size_map_for_current)
p2 = SpatialFeaturePlot(new.seurat.object, slot = "counts", features = "clusters", alpha = 0)
print(p1+p2)
dev.off()



