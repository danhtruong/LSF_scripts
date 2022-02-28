library(Seurat)
library(celda)
library(ggplot2)

#load data
MTL.integrated <- readRDS('inputs/2020-06-18 MTLv4 integrated harmony.RDS') 

#import counts from SCT assay
counts <- GetAssayData(MTL.integrated, assay = 'SCT', slot = 'counts')
counts <- as.matrix(counts)

cgs <- celdaGridSearch(counts,
    paramsTest = list(K = seq(10, 12), L = seq(18, 22)),
    cores = 16,
    model = "celda_CG",
    nchains = 2,
    maxIter = 100,
    verbose = TRUE,
    bestOnly = TRUE)

saveRDS(cgs, 'outputs/cgs.rds')

cgs <- resamplePerplexity(counts = counts,
    celdaList = cgs, resample = 5)

saveRDS(cgs, 'outputs/cgs_perplexity.rds')

#print pdf
saveRDS(cellSplit, 'outputs/MTLv4 gridsearch.rds')
pdf('outputs/cellSplit.pdf', width = 6, height = 4)
plotGridSearchPerplexity(celdaList = cgs) + theme(axis.text.x = element_text(size = 5, angle = 45))
dev.off()

