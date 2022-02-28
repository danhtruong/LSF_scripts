library(GSVA)
library(msigdbr)
library(Seurat)
library(dplyr)
setwd('/rsrch3/home/sarc_med_onco-rsch/dtruong4')

#load data
LPS_data.cell <- readRDS('inputs/2021-05-03 LPS_data.cell_integrated.rds')
msigdbr_df <- msigdbr(species = "Homo sapiens")
expr <- GetAssayData(LPS_data.cell, assay = 'SCT', slot ='data') %>% as.matrix()

#define gene sets
gene_sets <- c('GOMF_CHROMATIN_BINDING', 'MUELLER_PLURINET', 'KEGG_CELL_CYCLE')

msigdbr_df_gene_sets <- msigdbr_df %>% 
  filter(gs_name %in% gene_sets)

msigdbr_list = split(x = msigdbr_df_gene_sets$gene_symbol, f = msigdbr_df_gene_sets$gs_name)

#run GSVA

LPS_gsva <- gsva(expr = expr, 
gset.idx.list = msigdbr_list,
 method = 'ssgsea',
 kcdf = 'Gaussian', 
 parallel.sz = 28)
 
 saveRDS(LPS_gsva, 'outputs/2022-02-02 LPS_gsva.rds')