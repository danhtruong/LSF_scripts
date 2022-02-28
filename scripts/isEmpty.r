library(Seurat)
library(DropletUtils)
#DSRCT_PDX_data <- readRDS('inputs/2020-12-21 DSRCT_PDX_data.subset.rds')
ES_PDX_data <- readRDS( 'inputs/2020-12-21 ES_PDX_data.subset.rds')
OS_PDX_data <- readRDS('inputs/2021-01-05 OS_PDX_data.subset.rds')
LPS_data.cell.list <- readRDS('inputs/2020-10-12 LPS_data.cell.list.rds')
LPS_data.nuc.list <- readRDS('inputs/2020-10-12 LPS_data.nuc.list.rds')

#emptydrop function
emptyDropsResults <- function(seurat_object){
  my.counts <- GetAssayData(seurat_object, slot = 'counts')
  br.out <- barcodeRanks(my.counts)
  e.out <- emptyDrops(my.counts, lower = metadata(br.out)$knee)
}

OS_PDX_data <- SplitObject(OS_PDX_data, split.by = 'orig.ident')
e.out.list <- lapply(OS_PDX_data, emptyDropsResults)
saveRDS(e.out.list, '2021-01-22 OS_PDX_data.e.out.rds')

e.out.list <- lapply(LPS_data.cell.list, emptyDropsResults)
saveRDS(e.out.list, '2021-01-22 LPS_data.cell.e.out.rds')

e.out.list <- lapply(LPS_data.nuc.list, emptyDropsResults)
saveRDS(e.out.list, '2021-01-22 LPS_data.nuc.e.out.rds')

ES_PDX_data <- SplitObject(ES_PDX_data, split.by = 'orig.ident')
e.out.list <- lapply(ES_PDX_data, emptyDropsResults)
saveRDS(e.out.list, '2021-01-22 ES_PDX_data.e.out.rds')