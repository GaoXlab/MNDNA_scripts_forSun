###
# rm(list=ls())
library(stringr)
library(Rmisc)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(scales)
library(tidyverse)
library(reshape2)
library(RColorBrewer)
library(rtracklayer)
library(tidyverse)
library(GenomicRanges)
library(corrplot)
source('./plot_function.r')
load('../Human/02.mndna_enriched_peak_calling/data/ref.RData') # genome_size, feature_1000kb_anno, feature_100kb_anno
load('../Human/02.mndna_enriched_peak_calling/data/HD_deep_10sample.RData')
load('../Human/02.mndna_enriched_peak_calling/result/HD_deep_10sample_anno.RData')

dir.create('FigS3')

## Figure S3, pearson correlation
a=merge(HD12_60m_1000k_MN[,c('feature','median')], CN_smooth_r1_1000k, by.x='feature', by.y='region')
a=a[a$broadPeak==1, c('median','V2','overlappedGeneNumDensity')]; colnames(a) = c('MN-DNA\nread density','Chrom. size','Gene density')

cor2genomic <- cor(a[,])
cols <- rev(colorRampPalette(brewer.pal(9, "RdBu"))(500))
pdf('FigS3/FigS3.pdf', width=4, height=4)
corrplot(cor2genomic, method = 'color', col=rev(COL2('RdBu', 50)), addCoef.col = 'black', cl.cex = 0.6, number.cex=0.7)
dev.off()

