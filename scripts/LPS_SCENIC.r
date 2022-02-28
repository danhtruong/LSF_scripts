
library(Seurat)
library(SCENIC)

setwd('/rsrch3/home/sarc_med_onco-rsch/dtruong4/SCENIC')

### Load data
input <- readRDS('../inputs/2021-05-03 LPS_data.cell_integrated.rds')
exprMat <- GetAssayData(input, assay = 'SCT', slot = 'data')
exprMat <- as.matrix(exprMat)
cellInfo <- input@meta.data


### Initialize settings
org <- "hgnc" # or hgnc, or dmel
dbDir <- "../SCENIC/cisTarget_databases" # RcisTarget databases location
myDatasetTitle <- "SCENIC LPS DATA" # choose a name for your analysis

dbs <- c('10kb' = 'hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather', 
         '500bp' = 'hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather')
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=28) 

#scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
saveRDS(scenicOptions, file="../SCENIC/int/scenicOptions.Rds") 

### Co-expression network
genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,
                           minCountsPerGene=3*.01*ncol(exprMat),
                           minSamples=ncol(exprMat)*.01)
exprMat_filtered <- exprMat[genesKept, ]
runCorrelation(exprMat_filtered, scenicOptions)
#exprMat_filtered_log <- log2(exprMat_filtered+1) 

runGenie3(exprMat_filtered_log, scenicOptions)

scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions)
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log)
saveRDS(scenicOptions, file='../SCENIC/int/scenicOptions.Rds") # To save status