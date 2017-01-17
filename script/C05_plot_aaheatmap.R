#R code
library(ggplot2)
library(reshape)
library(scales)

plotheatmap <- function(y, graphname, h, w){
  colnames(y) <- c('pos','aa','value')
  p <- ggplot(data=y, aes(x=pos, y=aa, fill=value)) +
              geom_tile() +
              scale_fill_gradientn(colours=c("white", "yellow", "red"),
                         values=rescale(c(0, 0.5, 1)),
                         limits=c(0,0.78),
                         guide="colorbar") +
              theme_classic() +
              theme(panel.border = element_rect(colour = "black", fill=NA, size=4),
                    text = element_text(size=20)) +
              scale_x_discrete(expand=c(0,0)) + 
              scale_y_discrete(expand=c(0,0))
  ggsave(graphname,p, height=h, width=w)
  print (paste('File is written:',graphname,sep=' '))
  }

t <- read.table('data/AAFreqTable.tsv',header=1)

plotheatmap(t[,c(1,2,3)],'graph/HM_Input.png',6,3)
plotheatmap(t[,c(1,2,4)],'graph/HM_H1R1.png',6,3)
plotheatmap(t[,c(1,2,5)],'graph/HM_H1R2.png',6,3)
plotheatmap(t[,c(1,2,6)],'graph/HM_H1R3.png',6,3)
plotheatmap(t[,c(1,2,7)],'graph/HM_H3R1.png',6,3)
plotheatmap(t[,c(1,2,8)],'graph/HM_H3R2.png',6,3)
plotheatmap(t[,c(1,2,9)],'graph/HM_H3R3.png',6,3)
plotheatmap(t[,c(1,2,10)],'graph/HM_H5R1.png',6,3)
plotheatmap(t[,c(1,2,11)],'graph/HM_H5R2.png',6,3)
plotheatmap(t[,c(1,2,12)],'graph/HM_H5R3.png',6,3)
