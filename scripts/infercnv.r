library(Seurat)
library(infercnv)
require(dplyr)

options(expressions = 100000)

LPS_data.cell.list <- readRDS('../inputs/2021-05-03 LPS_data.cell_integrated.rds')
 gencode_v19_gene_pos <- read.delim("../inputs/gencode_v19_gene_pos.txt", header=FALSE, row.names=1)
 
LPS_data.cell.list <- SplitObject(LPS_data.cell.list, split.by = 'lab_id')
 
name <- names(LPS_data.cell.list)
lps_cnv <- list()
for (i in 1:length(LPS_data.cell.list))
{
  metadata <- LPS_data.cell.list[[i]]@meta.data[,'Annotations',drop=FALSE]
  ref <- unique(LPS_data.cell.list[[i]]@meta.data$Annotations)
  ref <- ref[!(ref == 'Fibroblasts' | ref == 'Adipocytes' | ref == 'Myocytes')]
  #ref <- c("Endothelial cells", 'Macrophages', 'NK cells', 'Monocytes', 'CD8+ T-cells', 'CD4+ T-cells')
  lps_cnv[[i]] = CreateInfercnvObject(raw_counts_matrix=LPS_data.cell.list[[i]]@assays[["SCT"]]@counts,
                                      annotations_file=metadata,
                                      gene_order_file=gencode_v19_gene_pos, ref_group_names=intersect(levels(metadata$Annotations), ref)) 
  
  
lps_cnv[[i]] = infercnv::run(lps_cnv[[i]],
                               cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                               out_dir=paste0('../outputs/CNV/', name[i]),
                               analysis_mode = 'subclusters',
                               cluster_by_groups=FALSE, 
                               denoise=TRUE,
                               HMM=TRUE,
							   HMM_type='i6',
							   no_plot = TRUE,
                               num_threads = 28)
}

#saveRDS(lps_cnv, '../outputs/2021-06-03 lps_cnv.rds')