#R code
library(ggplot2)
library(reshape)
library(scales)
require(cowplot)

FreqTable <- read.table('data/VariantFreqTable.tsv',header=1)
row.names(FreqTable) <- FreqTable$Variant
WT    <- FreqTable[which(FreqTable$Variant=='WT'),2:length(colnames(FreqTable))]
WT    <- WT[,which(!grepl('R4',colnames(WT)))]
WTfrac <- WT/1000000
Mutfrac <- 1-WTfrac
FracTable <- cbind(as.numeric(WTfrac),as.numeric(Mutfrac))
rownames(FracTable) <- colnames(WTfrac)
colnames(FracTable) <- c('WT','Mut')
FracTable <- melt(FracTable)
FracTable$X1 <- factor(FracTable$X1,levels=colnames(WT))
p <- ggplot(FracTable,aes(x=X1,y=value,fill=X2))+geom_bar(stat="identity")
ggsave('graph/WTFrac.png', p, height=3, width=4.5)
