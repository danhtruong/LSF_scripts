library(Seurat)
library(infercnv)
require(dplyr)

ES_patient <- readRDS('../inputs/ES_patient_data.rds')

gencode_v19_gene_pos <- read.delim("../inputs/gencode_v19_gene_pos.txt", header=FALSE, row.names=1)
 
metadata <- ES_patient@meta.data[,'Annotations', drop = FALSE]
ref <- c("Skeletal muscle" ,"Adipocytes", "Endothelial cells" , "CD4+ T-cells"  ) #Add reference 


cnv = CreateInfercnvObject(raw_counts_matrix=ES_patient@assays$RNA@counts,
                                      annotations_file=metadata,
                                      gene_order_file=gencode_v19_gene_pos, ref_group_names=intersect(unique(metadata$Annotations), ref)) 
  
  
cnv = infercnv::run(cnv,
                               cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                               out_dir='../outputs/ES/',
                               #analysis_mode = 'subclusters',
                               cluster_by_groups=TRUE, 
                               denoise=TRUE,
                               HMM=TRUE,
                               num_threads = 24)

saveRDS(cnv, '../outputs/2021-04-20 ES_cnv.rds')