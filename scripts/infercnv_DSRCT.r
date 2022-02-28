library(Seurat)
library(infercnv)
require(dplyr)
cell.list <- readRDS('../inputs/2021-04-30 DSRCT_patient_data.rds')
 gencode_v19_gene_pos <- read.delim("../inputs/gencode_v19_gene_pos.txt", header=FALSE, row.names=1)
 
cell.list <- SplitObject(cell.list, split.by = 'lab_id')
name <- names(cell.list)
cnv <- list()
for (i in 3:length(cell.list)
)
{
  metadata <- cell.list[[i]]@meta.data[,'Annotations',drop=FALSE]
  #ref <- unique(cell.list[[i]]@meta.data$Annotations)
  #ref <- ref[!(ref == 'Epith' | ref == 'Adipocytes' | ref == 'Chondrocytes')]
  ref <- c('CD4+ T-cells', 'Adipocytes')
  #ghetto fix for annotations, because some groups are below 1
  if (i == 2){
	levels = levels(metadata$Annotations)
	levels[levels == 'Smooth muscle'] = 'Epithelial cells'
	levels(metadata$Annotations) = levels
  }
  if (i == 3){
	levels = levels(metadata$Annotations)
	levels[levels == 'Epithelial cells'] = 'Smooth muscle'
	levels(metadata$Annotations) = levels
  }
  cnv[[i]] = CreateInfercnvObject(raw_counts_matrix=cell.list[[i]]@assays[["RNA"]]@counts,
                                      annotations_file=metadata,
                                      gene_order_file=gencode_v19_gene_pos, ref_group_names=intersect(unique(metadata$Annotations), ref)) 
  
  
cnv[[i]] = infercnv::run(cnv[[i]],
                               cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                               out_dir=paste0('../outputs/', name[i]),
                               #analysis_mode = 'subclusters',
                               cluster_by_groups=TRUE, 
                               denoise=TRUE,
                               HMM=TRUE,
                               num_threads = 24)
}

#saveRDS(cnv, '../outputs/2021-01-21 lps_cnv.rds')