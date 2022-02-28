library(Seurat)
library(celda)
library(ggplot2)

#load data
#input <- readRDS('inputs/2020-06-25 MTLv4 integrated harmony.RDS') 
input <- readRDS('inputs/2020-07-17 MTLv4 MSC.rds') 
input <- input[,!(input$orig.ident %in% 'Adipo 0H')]
#pr_graph_test_res <- readRDS('inputs/2020-06-27 pr_graph_test_res.rds')
#pr_deg_ids <- row.names(subset(pr_graph_test_res, q_value < 1e-3))


#import counts from SCT assay
counts <- GetAssayData(input, assay = 'SCT', slot = 'counts')
counts <- as.matrix(counts)#[pr_deg_ids,]

#find modules
celdaModel <- celda_CG(counts = counts, K = 4, L = 12, verbose = TRUE)
saveRDS(celdaModel, 'outputs/MTLv4_MSC_celdaModel_12_.rds')

#generate a heatmap
pdf('outputs/celdaHeatmap_MSC_12.pdf', width = 6, height = 4)
celdaHeatmap(counts = counts, celdaMod = celdaModel, nfeatures = 10)
dev.off()