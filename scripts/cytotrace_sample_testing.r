library(reticulate)
use_condaenv(condaenv = '/rsrch3/home/sarc_med_onco-rsch/dtruong4/.conda/envs/cytotrace', require = T)
py_config()

setwd('/rsrch3/home/sarc_med_onco-rsch/dtruong4')

library(Seurat)
library(dplyr)
LPS_data.cell <- readRDS('inputs/2021-05-03 LPS_data.cell_integrated.rds')

library(CytoTRACE)
LPS_data.cell.int <- GetAssayData(LPS_data.cell, assay = 'SCT', slot = 'count') %>% as.matrix()

results <- CytoTRACE(LPS_data.cell.int, ncores = 16, enableFast = T)

saveRDS(results, 'outputs/2022-02-02 LPS cytoTRACE.rds')
