
library(stringr)
library(ggplot2)
library(ggalt)
library(readxl)
library(stringr)
library(pheatmap)
library(stringr)
library(GenomicRanges)
source('./plot_function.r')


## crc features, overlapped with MN-DNA enriched regions, related to Fig2d
load('../Human/03.features_annotation/feature/features.RData')
crc_gr = GRanges(CRC0103)
mndna_gr = GRanges(MNDNApos6566)
overlapped_regions = length(unique(queryHits(findOverlaps(crc_gr, mndna_gr))))

df <- data.frame(
  name = c('1.overlapped','2.no_overlapped'), level = c('1','1'),
  values = c(round(100*(overlapped_regions/nrow(CRC0103)),2),round(100-100*(overlapped_regions/nrow(CRC0103)),2))
)
Fig2d <- ggplot(df, aes(x = level, y = values, fill = name, alpha = level)) +
        geom_col(width = 1, color = 'gray90', size = 1, position = position_stack(0.6)) + coord_polar(theta = 'y') +
        scale_alpha_manual(values = c('0' = 0.2, '1' = 0.8, '2' = 0.5), guide = F) +
        geom_text(aes(label = paste0(values, "%")), position = position_stack(vjust=0.5)) + 
        scale_x_discrete(breaks = NULL) + scale_y_continuous(breaks = NULL) + labs(x = NULL, y = NULL) + 
        scale_fill_manual(values=c("#1B7837","#D9D9D9"), na.translate = F) + theme_minimal()+theme(legend.position="none")
ggsave('./Fig2/Fig2d.pdf', Fig2d, width=3, height=3)



## examples, related to Fig2c
load('../Human/04.mndna_model/modelData/clinicalinfo.RData') # "cohort1_clinical" "cohort_INDP"
load('../Human/04.mndna_model/modelData/samplelist.RData') # "cohort1" "cohort_INDP2" "cohort_INDP3"
load('./Fig2/nor_100kb.RData')

plot_HDs = c("GLRHD0015","GLRHD0020","GLRHD0109","GLRHD0148","GLRHD0162","GLRHD0185","GLRHD0207","GLRHD0212","GLRHD0216","GLRHD0234")
plot_CRCs = c("GLRCRC0015","GLRCRC0047","GLRCRC0104","GLRCRC0133","GLRCRC0137","GLRCRC0142","GLRCRC0146","GLRCRC0162","GLRCRC0165","GLRCRC0167")

HD_medianInGroup <- MNdna_profiles_df1(nor_100kb, 'HD_med', intersect(colnames(nor_100kb), plot_HDs))
HD_medianInGroup <- as.data.frame(cbind(nor_100kb[, c('chr','start','end')], HD_medianInGroup)); HD_medianInGroup$median <- as.numeric(HD_medianInGroup$median)
CRC_medianInGroup <- MNdna_profiles_df1(nor_100kb, 'CRC_med', intersect(colnames(nor_100kb), plot_CRCs))
CRC_medianInGroup <- as.data.frame(cbind(nor_100kb[, c('chr','start','end')], CRC_medianInGroup)); CRC_medianInGroup$median <- as.numeric(CRC_medianInGroup$median)

# Feature1(6:162190000−162340000)
HD_medianInGroup_tmp = HD_medianInGroup[which((HD_medianInGroup$chr=='chr6')&(HD_medianInGroup$start>=161100000)&(HD_medianInGroup$start<=163300000)), ]
CRC_medianInGroup_tmp = CRC_medianInGroup[which((CRC_medianInGroup$chr=='chr6')&(CRC_medianInGroup$start>=161100000)&(CRC_medianInGroup$start<=163300000)), ]
start = 162190000
end = 162340000
p1 <- ggplot()+ 
    geom_vline(aes(xintercept=c(start, end)), colour="#FDC173", size=0.1, linetype='dashed') + 
    geom_rect(aes(xmin = start, xmax = end, ymin = 245, ymax = 456),fill = '#FDC173', alpha = 0.1)+
    geom_line(data=HD_medianInGroup_tmp, aes(x=start, y=median), size=0.5, color="#374E55FF")+
    geom_line(data=CRC_medianInGroup_tmp, aes(x=start, y=median), size=0.5, color="#E64B35FF")+
    xlab(str_c(as.character(chr),':',as.character(start),'-',as.character(end))) + ylab('Normalized read counts') +
    theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
ggsave('./Fig2/CRC0103.chr6:161100000-163300000.pdf',p1,width=6, height=3)

### Feature2(11:21130000−21220000)
HD_medianInGroup_tmp = HD_medianInGroup[which((HD_medianInGroup$chr=='chr11')&(HD_medianInGroup$start>=20600000)&(HD_medianInGroup$start<=22300000)), ]
CRC_medianInGroup_tmp = CRC_medianInGroup[which((CRC_medianInGroup$chr=='chr11')&(CRC_medianInGroup$start>=20600000)&(CRC_medianInGroup$start<=22300000)), ]
start = 21130000
end = 21220000
p1 <- ggplot()+ 
    geom_vline(aes(xintercept=c(start, end)), colour="#FDC173", size=0.1, linetype='dashed') + 
    geom_rect(aes(xmin = start, xmax = end, ymin = 218, ymax = 376),fill = '#FDC173', alpha = 0.1)+
    geom_line(data=HD_medianInGroup_tmp, aes(x=start, y=median), size=0.5, color="#374E55FF")+
    geom_line(data=CRC_medianInGroup_tmp, aes(x=start, y=median), size=0.5, color="#E64B35FF")+
    xlab(str_c(as.character(chr),':',as.character(start),'-',as.character(end))) + ylab('Normalized read counts') +
    theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
ggsave('./Fig2/CRC0103.chr11:20600000-22300000.pdf',p1,width=6, height=3)

### Feature3(2:213090000-213160000)
HD_medianInGroup_tmp = HD_medianInGroup[which((HD_medianInGroup$chr=='chr2')&(HD_medianInGroup$start>=212600000)&(HD_medianInGroup$start<=213500000)), ]
CRC_medianInGroup_tmp = CRC_medianInGroup[which((CRC_medianInGroup$chr=='chr2')&(CRC_medianInGroup$start>=212600000)&(CRC_medianInGroup$start<=213500000)), ]
start = 213090000
end = 213160000
p1 <- ggplot()+ 
    geom_vline(aes(xintercept=c(start, end)), colour="#FDC173", size=0.1, linetype='dashed') + 
    geom_rect(aes(xmin = start, xmax = end, ymin = 253, ymax = 518),fill = '#FDC173', alpha = 0.1)+
    geom_line(data=HD_medianInGroup_tmp, aes(x=start, y=median), size=0.5, color="#374E55FF")+
    geom_line(data=CRC_medianInGroup_tmp, aes(x=start, y=median), size=0.5, color="#E64B35FF")+
    xlab(str_c(as.character(chr),':',as.character(start),'-',as.character(end))) + ylab('Normalized read counts') +
    theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
ggsave('./Fig2/CRC0103.chr2:212600000-213500000.pdf',p1,width=6, height=3)

### Feature4(10:104090000−104250000)
HD_medianInGroup_tmp = HD_medianInGroup[which((HD_medianInGroup$chr=='chr10')&(HD_medianInGroup$start>=103900000)&(HD_medianInGroup$start<=104500000)), ]
CRC_medianInGroup_tmp = CRC_medianInGroup[which((CRC_medianInGroup$chr=='chr10')&(CRC_medianInGroup$start>=103900000)&(CRC_medianInGroup$start<=104500000)), ]
start = 104090000
end = 104250000
p1 <- ggplot()+ 
    geom_vline(aes(xintercept=c(start, end)), colour="#FDC173", size=0.1, linetype='dashed') + 
    geom_rect(aes(xmin = start, xmax = end, ymin = 348, ymax = 597),fill = '#FDC173', alpha = 0.1)+
    geom_line(data=HD_medianInGroup_tmp, aes(x=start, y=median), size=0.5, color="#374E55FF")+
    geom_line(data=CRC_medianInGroup_tmp, aes(x=start, y=median), size=0.5, color="#E64B35FF")+
    xlab(str_c(as.character(chr),':',as.character(start),'-',as.character(end))) + ylab('Normalized read counts') +
    theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
ggsave('./Fig2/CRC0103.chr10:103900000-104500000.pdf',p1,width=6, height=3)

pdf('./Fig2/Fig2.gene_locus.pdf')
getGeneLocus("chr6:161100000-163300000")
getGeneLocus("chr11:20600000-22300000")
getGeneLocus("chr2:212600000-213500000")
getGeneLocus("chr10:103900000-104500000")
dev.off()