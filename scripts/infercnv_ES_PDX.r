library(Seurat)
library(infercnv)
require(dplyr)

ES_PDX_data_msc <- readRDS('../inputs/ES_PDX_data_msc.rds')
ES_PDX_data_msc <- subset(ES_PDX_data_msc, downsample = 375)
 gencode_v19_gene_pos <- read.delim("../inputs/gencode_v19_gene_pos.txt", header=FALSE, row.names=1)


 
metadata <- ES_PDX_data_msc@meta.data[,'orig.ident', drop = FALSE]
ref <- c("MSC-8004L" ,"MSC-8011L")
ES_PDX_cnv_cnv = CreateInfercnvObject(raw_counts_matrix=ES_PDX_data_msc@assays$RNA@counts,
                                      annotations_file=metadata,
                                      gene_order_file=gencode_v19_gene_pos, ref_group_names=intersect(unique(metadata$orig.ident), ref)) 
  
  
ES_PDX_cnv_cnv = infercnv::run(ES_PDX_cnv_cnv,
                               cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                               out_dir='../outputs/ES-PDX/',
                               #analysis_mode = 'subclusters',
                               cluster_by_groups=TRUE, 
                               denoise=TRUE,
                               HMM=TRUE,
                               num_threads = 24)

saveRDS(ES_PDX_cnv_cnv, '../outputs/2021-10-05 ES_PDX_cnv.rds')