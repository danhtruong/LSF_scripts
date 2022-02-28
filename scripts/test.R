library(ggplot2)
pdf('moduleSplit.pdf', width = 6, height = 4)
p<- plotGridSearchPerplexity(celdaList = moduleSplit) + theme(axis.text.x = element_text(size = 5, angle= 45))
p
dev.off()