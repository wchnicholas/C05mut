#R code
library(ggplot2)
library(reshape)
library(scales)
require(cowplot)

t <- read.table('data/VariantFreqTable.tsv',header=1)
t <- t[which(t$Variant!='WT'),]
H3specialists <- t[which(log10(t$H1R3/sum(t$H1R3)) < -5 & log10(t$H3R3/sum(t$H3R3)) > -2.5),1]
H1specialists <- t[which(log10(t$H1R3/sum(t$H1R3)) > -2.4 & log10(t$H3R3/sum(t$H3R3)) < -4.5),1]
print('H3specialists:')
print(cat(as.character(H3specialists),sep="\n"))
print('H1specialists:')
print(cat(as.character(H1specialists),sep="\n"))
y <- data.frame(t$Variant,t$H1R3,t$H3R3)
colnames(y) <- c('Variant','H1R3','H3R3')
y$H1R3 <- y$H1R3/sum(y$H1R3)
y$H3R3 <- y$H3R3/sum(y$H3R3)
y[which(y$H1R3==0),]$H1R3 <- sort(y$H1R3[which(y$H1R3!=0)])[1]
y[which(y$H3R3==0),]$H3R3 <- sort(y$H3R3[which(y$H3R3!=0)])[1]
y$H1R3 <- log10(y$H1R3)
y$H3R3 <- log10(y$H3R3)
y_test <- y[sample(nrow(y), 10000), ]
p <- ggplot(y,aes(x=H1R3,y=H3R3)) +
       geom_rect(data=NULL,aes(xmin=-2.4,xmax=Inf,ymin=-Inf,ymax=-4.5), fill="#CC6666", alpha=0.02) +
       geom_rect(data=NULL,aes(xmin=-Inf,xmax=-4.75,ymin=-2.5,ymax=Inf), fill='#9999CC', alpha=0.02) +
       geom_point(size=0.5) + 
       xlim(-5.1,-0.9) +
       ylim(-5,-0.9)
ggsave('graph/H1H3specialists.png', p, height=4, width=4)
