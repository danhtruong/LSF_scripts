library(AUCell)
pred.matrix <- readRDS('pred.matrix.rds')
aucellRankings <- AUCell_buildRankings(pred.matrix,plotStats = FALSE, nCores = 28)
saveRDS(aucellRankings, 'aucellRankings.rds')