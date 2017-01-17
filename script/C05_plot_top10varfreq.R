#R code
library(ggplot2)
library(scales)
library(RColorBrewer)
require(cowplot)
library(readr)
library(dplyr)
library(tidyr)
library(reshape)
#library(extrafont)
#font_import(pattern="[A/a]rial")

CloneFreqPlot <- function (FreqTable,Pops,Top10Clones,h,w,graphname){
  pastel <- c(brewer.pal(12,"Set3"))
  pastel[2] <- '#CCCC00'
  rownames <- FreqTable$Variant
  FreqTable.melt <- FreqTable %>%
    gather(sample_name, frequency, -Variant)  %>%
    filter(Variant %in% Top10Clones | Variant == 'WT') %>%
    mutate(sample_name = factor(sample_name,levels=Pops))
  
  WT <- FreqTable.melt %>%
        filter(Variant == 'WT') %>%
	select(sample_name, frequency) %>%
	setNames(c('sample_name','WT_frequency')) 

  Top10Freq <- FreqTable.melt %>%
	inner_join(WT) %>%
        mutate(normalized_frequency = frequency / (1000000 - WT_frequency)) %>%
        filter(sample_name %in% Pops) %>%
        filter(Variant != 'WT') %>%
        arrange(normalized_frequency) %>%
        mutate(Variant=factor(Variant,levels=Top10Clones))

  p <- ggplot(data=Top10Freq, aes(x=sample_name, y=normalized_frequency, group=Variant, colour=Variant)) +
        geom_line() +
        scale_colour_manual(values=pastel) +
        ylim(0,0.15)
        theme(text=element_text(size=16, family="Arial"))

  ggsave(graphname, p, height=h, width=w)
  message('saved ', graphname)

  print (Pops[1])
  options(dplyr.print_max = Inf)
  print (filter(Top10Freq, grepl('R3',sample_name)),n=40,width = Inf)
  }

FreqTable <- read_tsv('data/VariantFreqTable.tsv')
graphheight <- 3
graphwidth  <- 6
H1R3Top10Clones <- c('VPGSGW','LPSSGW','VPGAGW','VASSGW','IPGSGW','VTGSGW','VATSGW','VVSSGW','VVGSGW','WPEIGF')
H3R3Top10Clones <- c('VPGAGW','WYVHLW','LPGGGW','YDPGGW','VPGSGW','VVSSGW','VVSAGW','VVDSGW','VASAGW','YEPAGW')
H5R3Top10Clones <- rev(c('LDAGDL','KCLWN_','DRPLAW','GHLHNW','ARELAY','NGCGRW','VVPEFW','RVLVRL','NPQEEL','YVHPQF'))
CloneFreqPlot(FreqTable, c('Input','H1R1','H1R2','H1R3'), H1R3Top10Clones, graphheight, graphwidth, 'graph/H1Top10Freq.png')
CloneFreqPlot(FreqTable, c('Input','H3R1','H3R2','H3R3'), H3R3Top10Clones, graphheight, graphwidth, 'graph/H3Top10Freq.png')
CloneFreqPlot(FreqTable, c('Input','H5R1','H5R2','H5R3'), H5R3Top10Clones, graphheight, graphwidth, 'graph/H5Top10Freq.png')
