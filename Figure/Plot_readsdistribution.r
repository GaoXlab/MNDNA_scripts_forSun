library(ggplot2)
library(ggpubr)
library(cowplot)
library(dplyr)
library(stringr)
library(reshape2)
source('./plot_function.r')
load('../Human/01.alignment_results/Summary_sequencing_reads.RData')
## plot fig1d/ fig2b/ sup fig2

## Figure 1d. gDNA reads distribution
df_deep = deep[, c('Sample', 'Intergenic%','Intron%','Exon%','Centromere%','mt%')]
df_deep = df_deep %>% mutate(sampleType = ifelse(str_detect(Sample, 'GLG'), 'gDNA', 'MNDNA'))
df_deep_melt = melt(df_deep, id.vars=c('Sample','sampleType'))
df3 <- data_summary(df_deep_melt, varname="value", groupnames=c("sampleType", "variable"))

p_intergenic <- ggplot(data=df3[df3$variable=='Intergenic%', ], aes(x=variable, y=value, fill=sampleType)) +ylim(0,100)+
    geom_bar(stat="identity", color="black", position=position_dodge())+ scale_fill_manual(values=c('#E14530', '#3E60AA'))+ 
    geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2, position=position_dodge(.9)) + mytheme3 
p_intron <- ggplot(data=df3[df3$variable=='Intron%', ], aes(x=variable, y=value, fill=sampleType)) +ylim(0,100)+
    geom_bar(stat="identity", color="black", position=position_dodge())+ scale_fill_manual(values=c('#E14530', '#3E60AA'))+ 
    geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2, position=position_dodge(.9))+ mytheme3
p_exon <- ggplot(data=df3[df3$variable=='Exon%', ], aes(x=variable, y=value, fill=sampleType)) +ylim(0,10)+
    geom_bar(stat="identity", color="black", position=position_dodge())+ scale_fill_manual(values=c('#E14530', '#3E60AA'))+ 
    geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2, position=position_dodge(.9))+ mytheme3
p_mt <- ggplot(data=df3[df3$variable=='mt%', ], aes(x=variable, y=value, fill=sampleType)) +ylim(0,0.5)+
    geom_bar(stat="identity", color="black", position=position_dodge())+ scale_fill_manual(values=c('#E14530', '#3E60AA'))+ 
    geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2, position=position_dodge(.9))+ mytheme3
p_centro <- ggplot(data=df3[df3$variable=='Centromere%', ], aes(x=variable, y=value, fill=sampleType)) +ylim(0,5)+
    geom_bar(stat="identity", color="black", position=position_dodge())+ scale_fill_manual(values=c('#E14530', '#3E60AA'))+ 
    geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2, position=position_dodge(.9))+ mytheme3
pg <- plot_grid(p_intergenic, p_intron, p_exon, p_centro, p_mt,ncol=5, align="hv",rel_widths=c(1.3,1,1,1))
ggsave('./Fig1/Fig1d.pdf', pg, height=2, width=4)


## Figure 2b. CRC cohort reads distribution
dir.create('Fig2')
cohort1_crc = cohort1[cohort1$label %in% c('HD','CRC'), ]
cohort1_crc$label = factor(cohort1_crc$label, levels=c('HD','CRC'))

colnames(cohort1_crc)[grep('Intergenic%', colnames(cohort1_crc))] = 'value_Intergenic'
p_intergenic <-  ggplot(data=cohort1_crc, aes(x=label, y=value_Intergenic, fill=label)) +
    geom_boxplot(outlier.size = 0.05, lwd=0.2)+scale_fill_manual(values=c("#374E55FF","#E64B35FF"))+ crctheme +
    stat_compare_means(label = "p.signif",method = "t.test",comparisons=list(c('HD','CRC')),label.x = 1.5,label.y = 54)+ylim(45,55)
colnames(cohort1_crc)[grep('Intron%', colnames(cohort1_crc))] = 'value_Intron'
p_intron <-  ggplot(data=cohort1_crc, aes(x=label, y=value_Intron, fill=label)) +
    geom_boxplot(outlier.size = 0.05, lwd=0.2)+scale_fill_manual(values=c("#374E55FF","#E64B35FF"))+ crctheme +
    stat_compare_means(label = "p.signif",method = "t.test",comparisons=list(c('HD','CRC')),label.x = 1.5,label.y = 49)+ylim(40,50)
colnames(cohort1_crc)[grep('Exon%', colnames(cohort1_crc))] = 'value_Exon'
p_exon <-  ggplot(data=cohort1_crc, aes(x=label, y=value_Exon, fill=label)) +
    geom_boxplot(outlier.size = 0.05, lwd=0.2)+scale_fill_manual(values=c("#374E55FF","#E64B35FF"))+ crctheme +
    stat_compare_means(label = "p.signif",method = "t.test",comparisons=list(c('HD','CRC')),label.x = 1.5,label.y = 2.7)+ylim(0,3)
colnames(cohort1_crc)[grep('Centromere%', colnames(cohort1_crc))] = 'value_Centromere'
p_centro <-  ggplot(data=cohort1_crc, aes(x=label, y=value_Centromere, fill=label)) +
    geom_boxplot(outlier.size = 0.05, lwd=0.2)+scale_fill_manual(values=c("#374E55FF","#E64B35FF"))+ crctheme +
    stat_compare_means(label = "p.signif",method = "t.test",comparisons=list(c('HD','CRC')),label.x = 1.5,label.y = 2.7)+ylim(0,3)
colnames(cohort1_crc)[grep('mt%', colnames(cohort1_crc))] = 'value_mt'
p_mt <-  ggplot(data=cohort1_crc, aes(x=label, y=value_mt, fill=label)) +
    geom_boxplot(outlier.size = 0.05, lwd=0.2)+scale_fill_manual(values=c("#374E55FF","#E64B35FF"))+ crctheme +
    stat_compare_means(label = "p.signif",method = "t.test",comparisons=list(c('HD','CRC')),label.x = 1.5,label.y = 0.09)+ylim(0,0.1)
pg <- plot_grid(p_intergenic, p_intron, p_exon,p_centro, p_mt, ncol=5, align="hv",rel_widths=c(1,1,1,1))#p_centro, ncol=4,
ggsave('Fig2/Fig2b.pdf', pg, width=6, height=2)


### Supplementary figures
dir.create('FigS2')
## sup fig2a
load('../Human/02.mndna_enriched_peak_calling/data/HD_deep_downsamp_10sample.RData')
individual_list = c("GLHD0001","GLHD0002","GLHD0003","GLHD0004","GLHD0005","GLHD0006","GLHD0007","GLHD0008","GLHD0009","GLHD0010")
df_merge = c()
d_120m = deep_downsamp[['120m']]
for(downsamp in c("20m","30m","40m","50m","60m","70m","80m","90m","100m","110m","120m")){
    d_tmp = deep_downsamp[[downsamp]]
    MNDNAcor = c()
    gDNAcor = c()
    MNDNA2gDNA = c()
    for(individual in individual_list){
        mnDNAid = gsub('^GL', 'GLR', individual)
        gDNAid = gsub('^GL', 'GLG', individual)
        mnDNA_fea = merge(d_120m[, c('feature', mnDNAid)], d_tmp[, c('feature', mnDNAid)], by='feature')
        gDNA_fea = merge(d_120m[, c('feature', gDNAid)], d_tmp[, c('feature', gDNAid)], by='feature')
        # spearman cor
        MNDNAcor = c(MNDNAcor, cor(mnDNA_fea[,2:3], method='spearman')[1,2])
        gDNAcor = c(gDNAcor, cor(gDNA_fea[,2:3], method='spearman')[1,2])
        MNDNA2gDNA = c(MNDNA2gDNA, cor(d_tmp[, mnDNAid], d_tmp[, gDNAid], method='spearman'))
    }
    MNDNAcor = as.data.frame(MNDNAcor); MNDNAcor$label = 'MNDNAcor'; colnames(MNDNAcor)[1] = 'value'
    gDNAcor = as.data.frame(gDNAcor); gDNAcor$label = 'gDNAcor'; colnames(gDNAcor)[1] = 'value'
    MNDNA2gDNA = as.data.frame(MNDNA2gDNA); MNDNA2gDNA$label = 'MNDNA2gDNA'; colnames(MNDNA2gDNA)[1] = 'value'
    df = as.data.frame(rbind(MNDNAcor, gDNAcor, MNDNA2gDNA))
    df$downsamp = downsamp
    df_merge = rbind(df_merge, df)
}
df1 = as.data.frame(df_merge); 
df2 <- data_summary(as.data.frame(df1), varname='value', groupnames=c('downsamp','label'))
df2$downsamp = factor(df2$downsamp, levels=rev(c("20m","30m","40m","50m","60m","70m","80m","90m","100m","110m","120m")))
df2$label = factor(df2$label, levels=c('MNDNAcor','gDNAcor','MNDNA2gDNA'))
figs2a <- ggplot(df2, aes(x= downsamp, y = value, group=label, color=label)) +
        geom_line(position = position_dodge(0.1), linetype="dashed") +
        geom_point(position = position_dodge(0.1), size=1) + geom_errorbar(aes(ymin = value - sd, ymax = value + sd), width = 0.05, position = position_dodge(0.1)) + 
        scale_color_manual(values=c("#E64B35FF","#3E60AA","#DF8F44FF")) +
        theme_bw()+theme(legend.position="bottom",legend.title=element_blank(),axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1)) + 
        scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.25)) + xlab('') + ylab('Correlation to initial MN-DNA\nread density at 6x coverage') + xlab('Subsampling reads (# million)')
# ggsave('./FigS2/FigS2a.pdf',p1,width=3.2,height=2.5)

## sup fig2d
load('../Human/02.mndna_enriched_peak_calling/data/HD_replicants.RData')
figs2b <- ggscatter(HD_replicants, x = 'GLRHD1003', y = 'GLRHD1004', size=1,
          add = "reg.line", conf.int = TRUE, fill = "#E64B35FF",color = "#E64B35FF",
          cor.coef = TRUE, cor.method = "spearman",add.params = list(color="#374E55FF", type="dashed"), 
          xlab = "MN-DNA (replicant 1)", ylab = "MN-DNA (replicant 2)")
# ggsave('./FigS2/FigS2d.pdf',p1,width=3.2,height=2.5)

## sup fig2c
library(RColorBrewer)
load('../Human/01.alignment_results/Exp_sampleinfo.RData')
load('../Human/01.alignment_results/Summary_sequencing_reads.RData')

exp = merge(exp, Exp_sampleinfo[,c('Sample','Sample.Type')], by='Sample')
exp$Sample.Type = gsub('MN-DNA\\ ', '', gsub(', technical replicate','',exp$Sample.Type))
exp$Sample.Type = gsub("gDNA\\ \\(\\+\\ 200pg\\ gDNA)",'gDNA',exp$Sample.Type)
gDNA_contamination <- c('+ 0pg gDNA','+ 0.6pg gDNA','+ 1.2pg gDNA','+ 1.8pg gDNA','+ 3pg gDNA','+ 6pg gDNA','+ 12pg gDNA','+ 18pg gDNA','+ 24pg gDNA','+ 30pg gDNA','+ 60pg gDNA','+ 90pg gDNA','+ 120pg gDNA','+ 150pg gDNA','+ 180pg gDNA','gDNA')
exp$Sample.Type <- factor(exp$Sample.Type, levels=gDNA_contamination)
colnames(exp)[grep('mt%', colnames(exp))] = 'value_mt'
figs2c <- ggplot(data = exp, aes(x = Sample.Type, y = value_mt , fill=Sample.Type)) + ##
        geom_boxplot(outlier.size = 0.003, lwd=0.1)+theme_bw()+ geom_jitter(width = 0.3, size = 0.003)+#
        geom_hline(yintercept = 0.04,linetype = "dashed", color='grey')+ylab('Proportion of MN-DNA\nor leukocyte DNA\nmapped to MT regions (%)')+ xlab('Added leukocyte DNA (pg)')+ ylim(0,0.15)+
        scale_fill_manual(values=c(rev(brewer.pal(7,"Reds")), brewer.pal(9,"Blues")))+ mytheme2
# ggsave('./FigS2/FigS2c.pdf',figs2c,width=3.2,height=2.5)

## sup fig2d
load('../Human/01.alignment_results/Summary_downsamp_sequencing_reads.RData') 
Summary_downsamp_sequencing_reads$mt_perc <- 100*Summary_downsamp_sequencing_reads$mt_perc
Summary_downsamp_sequencing_reads$label = factor(Summary_downsamp_sequencing_reads$label, levels=c('MN DNA','leukocyte DNA'))
figs2d <- ggplot(data = Summary_downsamp_sequencing_reads[(Summary_downsamp_sequencing_reads$downsampling=='20m'), ], aes(x = label, y = mt_perc )) + ##
        geom_boxplot(outlier.colour = NA, fill=c('#E64B35FF', '#3E60AA'), lwd=0.1)+geom_jitter(width = 0.3, size = 0.003)+
        geom_hline(yintercept = 0.04,linetype = "dashed", color='grey')+
        theme_classic()+ theme(legend.position='none',axis.title.x=element_blank(), axis.text.x = element_text(angle = 45,hjust = 1))+ylim(0,0.15)+ylab('Proportion of MN-DNA\nor leukocyte DNA\nmapped to MT regions (%)')
# ggsave('./FigS2/FigS2d.pdf',g,width=2.5,height=1.9) ###sup用

## sup fig2f
map2mouse <- c(0.9606,0.9257,0.9767,0.9608,0.0003,0.0002)
map2hus <- c(0.0032,0.0046,0.0009,0.0012,0.9980,0.9981 )
label <- c('MN-DNA','MN-DNA','gDNA','gDNA','H460 cell line','H460 cell line')
df = as.data.frame(cbind(map2mouse, map2hus, label))
df$map2mouse = as.numeric(df$map2mouse)
df$map2hus = as.numeric(df$map2hus)
df$label = factor(df$label, levels=c('MN-DNA','gDNA','H460 cell line'))
figs2f <-ggplot(df, aes(x=map2mouse, y=map2hus, color=label)) + geom_point(shape=18)+ylim(0,1)+xlim(0,1)+
        geom_hline(yintercept=0, linetype="dashed", color='grey', lwd=0.5)+ geom_vline(xintercept=0, linetype="dashed", color='grey', lwd=0.5)+
        xlab('reads mapped to\nmouse genome (mm10)')+ylab('reads mapped to\nhuman genome (GRCh38)')+
        scale_color_manual(values=c("#E64B35FF","#3E60AA","#DF8F44FF"))+theme_classic()
# ggsave('./FigS2/FigS2f.pdf',figs2f,width=3,height=2)



blank_plot <- ggplot() + theme_void()
g <- plot_grid(plot_grid(figs2a,figs2b,figs2c, nrow=1),
               plot_grid(figs2d,figs2f,blank_plot, nrow=1, rel_widths=c(0.5, 1)), nrow=2, rel_heights=c(1.2,1))
ggsave('./FigS2/FigS2.pdf',g,width=10,height=6)


