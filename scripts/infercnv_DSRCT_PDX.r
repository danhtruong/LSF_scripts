library(Seurat)
library(infercnv)
require(dplyr)

DSRCT_PDX_data_msc <- readRDS('../inputs/DSRCT_PDX_data_msc.rds')
DSRCT_PDX_data_msc <- subset(DSRCT_PDX_data_msc, downsample = 375)
 gencode_v19_gene_pos <- read.delim("../inputs/gencode_v19_gene_pos.txt", header=FALSE, row.names=1)


 
metadata <- DSRCT_PDX_data_msc@meta.data[,'orig.ident', drop = FALSE]
ref <- c("MSC-8004L" ,"MSC-8011L")
DSRCT_PDX_cnv = CreateInfercnvObject(raw_counts_matrix=DSRCT_PDX_data_msc@assays$RNA@counts,
                                      annotations_file=metadata,
                                      gene_order_file=gencode_v19_gene_pos, ref_group_names=intersect(unique(metadata$orig.ident), ref)) 
  
  
DSRCT_PDX_cnv = infercnv::run(DSRCT_PDX_cnv,
                               cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                               out_dir='../outputs/DSRCT-PDX/',
                               #analysis_mode = 'subclusters',
                               cluster_by_groups=TRUE, 
                               denoise=TRUE,
                               HMM=TRUE,
                               num_threads = 24)

saveRDS(DSRCT_PDX_cnv, '../outputs/2021-10-05 DSRCT_PDX_cnv.rds')