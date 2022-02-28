library(Seurat)
library(celda)
library(ggplot2)

#load data
#input <- readRDS('inputs/2020-06-25 MTLv4 integrated harmony.RDS') 
input <- readRDS('inputs/2020-07-17 MTLv4 MSC.rds') 
moduleSplit <- readRDS('outputs/MTLv4 MSC 25 moduleSplit.rds')

#import counts from SCT assay
counts <- GetAssayData(input, assay = 'SCT', slot = 'counts')
counts <- as.matrix(counts)

moduleSplitSelect <- subsetCeldaList(moduleSplit, params = list(L = 12))
cellSplit <- recursiveSplitCell(counts = counts,
    initialK = 1,
    maxK = 5,
    yInit = clusters(moduleSplitSelect)$y)

#print pdf
saveRDS(cellSplit, 'outputs/MTLv4 MSC cellSplit.rds')
pdf('outputs/cellSplit.pdf', width = 6, height = 4)
plotGridSearchPerplexity(celdaList = cellSplit) + theme(axis.text.x = element_text(size = 5, angle = 45))
dev.off()

