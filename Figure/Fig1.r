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
library(reshape2)
source('./plot_function.r')
load('../Human/02.mndna_enriched_peak_calling/data/ref.RData') # genome_size, feature_1000kb_anno, feature_100kb_anno
load('../Human/02.mndna_enriched_peak_calling/data/HD_deep_10sample.RData') # "HD10_60m_1000k_MN"   "HD10_60m_100k_MN"    "HD10_60m_1000k_gDNA" "HD10_60m_100k_gDNA"
load('../Human/02.mndna_enriched_peak_calling/result/HD_deep_10sample_anno.RData') # "CN_smooth_r1_1000k" "CN_smooth_r1_100k"

dir.create('Fig1')

## Fig1e: genomic distribution across the genome
p0 <- ggplot() + #
    geom_hline(yintercept=1, size=0.3, linetype='dashed', color='grey')+#, colour="#67000D", size=0.1, linetype='dashed') + 
    geom_ribbon(data=HD10_60m_1000k_MN, aes(x=start, ymin = min, ymax = max), alpha=0.6, fill="#E64B35FF")+
    geom_line(data=HD10_60m_1000k_MN, aes(x=start, y=median, group=1), color="#E64B35FF", size=0.3)+
    theme_bw()+facet_grid(.~arm, switch="x", space="free", scales="free")+mytheme+ylim(0.5,2)+theme(strip.text.x = element_blank())+ylab('Median read counts')
p2 <- ggplot() + #
    geom_hline(yintercept=1, size=0.3, linetype='dashed', color='grey')+#, colour="#67000D", size=0.3, linetype='dashed') + 
    geom_ribbon(data=HD10_60m_1000k_gDNA, aes(x=start, ymin = min, ymax = max), alpha=0.6, fill="#3E60AA")+
    geom_line(data=HD10_60m_1000k_gDNA, aes(x=start, y=median, group=1), color="#3E60AA", size=0.1)+
    theme_bw()+facet_grid(.~arm, switch="x", space="free", scales="free")+mytheme+ylim(0.5,2)+theme(strip.text.x = element_blank())+ylab('Median read counts')
p1_anno <- ggplot(CN_smooth_r1_1000k)+
        geom_tile(aes(x=start, y=broadPeak_y, fill = broadPeak) ,colour = NA)+
        scale_fill_gradient2(low = "white", high = "red4")+mytheme+facet_grid(.~arm,switch="x", space="free_x", scales="free")+ylab('enriched')
g <- plot_grid(p0, p2, p1_anno, ncol=1,rel_heights=c(1,1,0.6))
ggsave('Fig1/Fig1e.pdf', g, width=8, height=2.3)


## Fig1f: chr1:198000001-200100001 zoom out
p0 <- ggplot() + #
    geom_hline(yintercept=1, size=0.3, linetype='dashed', color='grey')+#, colour="#67000D", size=0.1, linetype='dashed') + 
    geom_ribbon(data=HD10_60m_100k_MN[(HD10_60m_100k_MN$chromosome==1)&(HD10_60m_100k_MN$start>198000001)&(HD10_60m_100k_MN$end<200100001), ], aes(x=start, ymin = min, ymax = max), alpha=0.5, fill="#E64B35FF")+
    geom_line(data=HD10_60m_100k_MN[(HD10_60m_100k_MN$chromosome==1)&(HD10_60m_100k_MN$start>198000001)&(HD10_60m_100k_MN$end<200100001), ], aes(x=start, y=median, group=1), color="#E64B35FF", size=0.3)+
    theme_bw()+facet_grid(.~arm, switch="x", space="free", scales="free")+mytheme+ylim(0,4)+theme(strip.text.x = element_blank())+ylab('Median read counts')
p2 <- ggplot() + #
    geom_hline(yintercept=1, size=0.3, linetype='dashed', color='grey')+#, colour="#67000D", size=0.3, linetype='dashed') + 
    geom_ribbon(data=HD10_60m_100k_gDNA[(HD10_60m_100k_gDNA$chromosome==1)&(HD10_60m_100k_gDNA$start>198000001)&(HD10_60m_100k_gDNA$end<200100001), ], aes(x=start, ymin = min, ymax = max), alpha=0.5, fill="#3E60AA")+
    geom_line(data=HD10_60m_100k_gDNA[(HD10_60m_100k_gDNA$chromosome==1)&(HD10_60m_100k_gDNA$start>198000001)&(HD10_60m_100k_gDNA$end<200100001), ], aes(x=start, y=median, group=1), color="#3E60AA", size=0.3)+
    theme_bw()+facet_grid(.~arm, switch="x", space="free", scales="free")+mytheme+ylim(0,4)+theme(strip.text.x = element_blank())+ylab('Median read counts')
p1_anno <- ggplot(CN_smooth_r1_100k[(CN_smooth_r1_100k$chr=='chr1')&(CN_smooth_r1_100k$start>198000001)&(CN_smooth_r1_100k$end<200100001), ])+
        geom_tile(aes(x=start, y=broadPeak_y, fill = broadPeak) ,colour = NA)+
        scale_fill_gradient2(low = "white", high = "red4")+mytheme+facet_grid(.~arm,switch="x", space="free_x", scales="free")+ylab('enriched')
g = plot_grid(p0, p2, p1_anno, ncol=1,rel_heights=c(1,1,0.6))
ggsave('Fig1/Fig1f_zoomout.pdf', g, width=4, height=2.3)


## Fig1g: 
order_row = c("weak.enhancers","transcription","enhancers","quescient","transcribed.and.enhancer","polycomb.repressed","acetylations","promoters","weak.transcription","bivalent.promoters","exon","HET","znf","TSS","others","DNase")
mndna_enriched = read.table('../Human/03.features_annotation/result/MNDNApos6566.grouped_output.txt', sep='\t', head=TRUE, row.names=1)

df=mndna_enriched[apply(mndna_enriched,1,sum)<101, ]
tmp = t(mndna_enriched)[intersect(order_row, rownames(t(mndna_enriched))), ]
mndna_enriched_summary = cal_count(tmp)

num = mndna_enriched_summary[[2]]
active = sum(num[(num$Var2=='weak.enhancers')|(num$Var2=='transcription')|(num$Var2=='enhancers')|(num$Var2=='transcribed.and.enhancer')|
      (num$Var2=='acetylations')|(num$Var2=='promoters')|(num$Var2=='weak.transcription')|(num$Var2=='bivalent.promoters')|
      (num$Var2=='exon')|(num$Var2=='TSS')|(num$Var2=='DNase'),'value'])
noactive = sum(num$value) - active
noanno = nrow(mndna_enriched) - active - noactive
total_perc = as.data.frame(cbind(c(active,noactive,noanno),c('active','no-active','no-anno')))
total_perc$V1 = as.numeric(total_perc$V1)
total_perc$perc = (total_perc$V1/nrow(mndna_enriched))*100
total_perc$label = 'perc'
total_perc$V2 = factor(total_perc$V2, levels=c('no-anno','no-active','active'))
fig1g <- ggplot(total_perc, aes(x=label,y=perc,fill=V2))+geom_bar(stat="identity",alpha=0.5)+
        scale_fill_manual(values=c('grey','#4D7FB9','red4'))+ylab('Proportion of overlapped>50%\nwithin annotated chromatin states')+
        theme_bw()+theme(legend.position='none', panel.grid.major=element_blank(), panel.grid.minor=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank())
# ggsave('Fig1/Fig1g.pdf',fig1h, width=2, height=3)



random = read.table('../Human/03.features_annotation/result/random/all10times.cal_perc.log', sep='\t', head=TRUE, row.names=1)

random_melt = melt(t(random))
random_melt$Var2 = factor(random_melt$Var2, levels=rev(order_row))
random_melt$label = 'random'    
random_melt$value_flip = -random_melt$value
df2 <- data_summary(random_melt[random_melt$label=='random',], varname="value_flip", groupnames=c("Var2", "label"))
df2$Var2=as.factor(df2$Var2)

p_random <- ggplot(df2, aes(x=Var2, y=value_flip), fill='grey', alpha=0.5) + 
    geom_bar(stat = "identity", position = "identity") + geom_errorbar(aes(ymin=value_flip-sd, ymax=value_flip+sd), width=.2, position=position_dodge(.9)) +
    coord_flip() + theme_bw() +theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.title.x=element_blank()) + xlab('Proportion of regions\nwith chromatin marks (%)')
p_mnDNA <- ggplot(mndna_enriched_summary[[1]], aes(x=Var2, y=value)) + 
    geom_bar(stat = "identity", position = "identity", fill='red4', alpha=0.5) +
    coord_flip() + theme_bw() + 
    theme(axis.text.y=element_blank(),axis.title.x=element_blank(), axis.ticks.y=element_blank(),panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + xlab('') + ylim(0,0.3)

fig1g_m = plot_grid(fig1g, p_random,p_mnDNA, ncol=3, rel_widths=c(0.5, 1.05, 0.55))
ggsave('Fig1/Fig1g.pdf', width=6.5, height=3)

