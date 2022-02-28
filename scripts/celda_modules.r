library(Seurat)
library(celda)
library(ggplot2)

#load data
#input <- readRDS('inputs/2020-06-25 MTLv4 integrated harmony.RDS') 
input <- readRDS('inputs/2020-07-17 MTLv4 MSC.rds') 
#pr_graph_test_res <- readRDS('inputs/2020-06-27 pr_graph_test_res.rds')
#pr_deg_ids <- row.names(subset(pr_graph_test_res, q_value < 1e-3))




#import counts from SCT assay
counts <- GetAssayData(input, assay = 'SCT', slot = 'counts')
counts <- as.matrix(counts)#[pr_deg_ids,]

#find optimal number of modules
moduleSplit <- recursiveSplitModule(counts = counts,initialL = 2, maxL = 25) #normally 50 but this is a small dataset

#print pdf
saveRDS(moduleSplit, 'outputs/MTLv4 MSC 25 moduleSplit.rds')
pdf('outputs/moduleSplit.pdf', width = 6, height = 4)
plotGridSearchPerplexity(celdaList = moduleSplit) + theme(axis.text.x = element_text(size = 5, angle = 45))
dev.off()
