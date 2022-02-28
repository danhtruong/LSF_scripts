library(cisTopic)
library(Seurat)
library(data.table)
MTL_ATAC.subset <-  readRDS('inputs/2020-07-20 MTLv4_ATAC_subset.rds')


mtx <- MTL_ATAC.subset@assays$merged_peaks@counts
rnames = data.table('region' = rownames(mtx))
tmp = tidyr::separate(rnames, col = 'region', into = c('chr', 'start', 'end'))
rnames = paste0(tmp$chr, ':', tmp$start, '-', tmp$end)
rownames(mtx) = rnames

cisTopicObject <- createcisTopicObject(mtx, project.name='scATAC')

cisTopicObject <- runWarpLDAModels(cisTopicObject, topic= c(10:25, 30, 35, 40, 45, 50), 
seed=987, nCores=28, iterations = 500, addModels=FALSE)
saveRDS(cisTopicObject, 'outputs/2020-07-18 cistopicobject.rds')