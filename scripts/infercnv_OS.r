library(Seurat)
library(infercnv)
require(dplyr)
OS_PDX_data_OB <- readRDS('inputs/OS_PDX_data_OB.rds')
OS_PDX_data_OB <- subset(OS_PDX_data_OB, downsample = 1000)
 gencode_v19_gene_pos <- read.delim("inputs/gencode_v19_gene_pos.txt", header=FALSE, row.names=1)
 
metadata <- OS_PDX_data_OB@meta.data[,'orig.ident', drop = FALSE]
ref <- c('Osteoblasts')
OS_cnv = CreateInfercnvObject(raw_counts_matrix=OS_PDX_data_OB@assays$RNA@counts,
                                      annotations_file=metadata,
                                      gene_order_file=gencode_v19_gene_pos, ref_group_names=intersect(unique(metadata$orig.ident), ref)) 
  
  
OS_cnv = infercnv::run(OS_cnv,
                               cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                               out_dir='outputs/OS/',
                               analysis_mode = 'subclusters',
                               cluster_by_groups=TRUE, 
                               denoise=TRUE,
                               HMM=TRUE,
                               num_threads = 24)

saveRDS(OS_cnv, 'outputs/2021-03-16 OS_cnv.rds')
