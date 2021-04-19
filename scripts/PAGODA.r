library(Seurat)
library(scde)
library(parallel)
require(dplyr)
LPS_data.cell.list <- readRDS('inputs/2020-10-12 LPS_data.cell.list.rds')
n.cores <- detectCores()
print(n.cores)
cd <- LPS_data.cell.list[[1]]@assays$SCT@counts %>% as.matrix()
cd <-apply(cd,2,function(x) {storage.mode(x) <- 'integer'; x})


knn <- knn.error.models(cd, k = ncol(cd)/30, n.cores = n.cores, min.count.threshold = 1, min.nonfailed = 50, max.model.plots = 10)

varinfo <- pagoda.varnorm(knn, counts = cd, trim = 3/ncol(cd), max.adj.var = 5, n.cores = n.cores, plot = TRUE)
varinfo <- pagoda.subtract.aspect(varinfo, colSums(cd[, rownames(knn)]>0))


library(org.Hs.eg.db)
# translate gene names to ids
ids <- unlist(lapply(mget(rownames(cd), org.Hs.egALIAS2EG, ifnotfound = NA), function(x) x[1]))
rids <- names(ids); names(rids) <- ids 
# convert GO lists from ids to gene names
gos.interest <- unique(c(ls(org.Hs.egGO2ALLEGS)[1:100],"GO:0022008","GO:0048699", "GO:0000280")) 
go.env <- lapply(mget(gos.interest, org.Hs.egGO2ALLEGS), function(x) as.character(na.omit(rids[x]))) 
go.env <- clean.gos(go.env) # remove GOs with too few or too many genes
go.env <- list2env(go.env) # convert to an environment

pwpca <- pagoda.pathway.wPCA(varinfo, go.env, n.components = 1, n.cores = n.cores)

clpca <- pagoda.gene.clusters(varinfo, trim = 7.1/ncol(varinfo$mat), n.clusters = 50, n.cores = n.cores, plot = TRUE)

save.image(file = "PAGODA.RData")