
library(Seurat)
library(LandSCENT)
library(SingleCellExperiment)

setwd('/rsrch3/home/sarc_med_onco-rsch/dtruong4')

#load data
input <- readRDS('inputs/2020-06-25 MTLv4 integrated harmony.RDS')
input.sce <- as.SingleCellExperiment(MTL.integrated)
pheno.v = input.sce$anno


#make sure the gene name matches  with PPI
require(AnnotationDbi)
require(org.Hs.eg.db)
example.m <- as.matrix(assay(input.sce, i = "logcounts"))
anno.v <- mapIds(org.Hs.eg.db, keys = rownames(example.m), keytype = "SYMBOL", 
                 column = "ENTREZID", multiVals = "first")
unique_anno.v <- unique(anno.v)
example_New.m <- matrix(0, nrow = length(unique_anno.v), ncol = dim(example.m)[2])
for (i in seq_len(length(unique_anno.v))) {
  tmp <- example.m[which(anno.v == unique_anno.v[i]) ,]
  if (!is.null(dim(tmp))) {
    tmp <- colSums(tmp) / dim(tmp)[1]
  }
  example_New.m[i ,] <- example_New.m[i ,] + tmp
}
rownames(example_New.m) <- unique_anno.v
colnames(example_New.m) <- colnames(example.m)
example_New.m <- example_New.m[-which(rownames(example_New.m) %in% NA) ,]
Example.m <- example_New.m

#remove all zeros
Example.m[Example.m ==0 ] <- 0.1

#load ppi network
data(net17Jan16.m)
Integration.l <- DoIntegPPI(exp.m = Example.m, ppiA.m = net17Jan16.m)
str(Integration.l)

#compute SR
SR.o <- CompSRana(Integration.l, local = TRUE, mc.cores = 16)
saveRDS(SR.o,'outputs/2020-07-06 MTL signaling_entropy.rds')
#SR.o <- readRDS('2020-06-12 GSE113253 signaling_entropy.rds')

#infer potency
InferPotency.o <- InferPotency(SR.o, pheno.v = pheno.v, maxPS =15)

#infer landmarks
InferLandmark.o <- InferLandmark(InferPotency.o, pheno.v = pheno.v,
                                 reduceMethod = "PCA", clusterMethod = "PAM",
                                 k_pam = 15)
saveRDS(InferLandmark.o, 'outputs/2020-07-06 MTL InferLandmark.o.RDS')		


LandSR.o <- Plot_LandSR(InferLandmark.o, coordinates = Embeddings(input, 'umap'), bty = "f", PDF = TRUE, lighting = TRUE, scale_z = TRUE)						
