library(pheatmap)
library(GenomicRanges)
library(stringr)
library(reshape)
library(ggplot2)
library(ggpubr)

dir.create('./Fig2/')
dir.create('./FigS6/')


## crc features, chromHMM annotation, related to FigS6
CRC_246 = read.table('../Human/03.features_annotation/result/CRC0103.grouped_output.txt', head=TRUE, sep='\t', row.names=1)
CRC_246_chromHMM = CRC_246[names(sort(apply(CRC_246,2,sum),  decreasing = TRUE))]
chr_order = c("quescient","polycomb.repressed","HET","transcription","enhancers","weak.enhancers","acetylations","weak.transcription","transcribed.and.enhancer","exon","promoters","others","znf","bivalent.promoters","TSS","DNase" )
CRC_246_chromHMM = CRC_246[chr_order]

pdf('./FigS6/FigS6.chromHMM.pheatmap.pdf', width=5, height=6)
pheatmap(CRC_246_chromHMM, cluster_rows=TRUE,cluster_cols=FALSE,scale='none',show_rownames=FALSE, #clustering_distance_rows='correlation',
         color = c(colorRampPalette(colors = c("white","#cb181d","firebrick3"))(250)))
dev.off()


## crc features, normalized read counts distribution, related to FigS6
load('../Human/04.mndna_model/modelData/clinicalinfo.RData') # "cohort1_clinical" "cohort_INDP"
load('../Human/04.mndna_model/modelData/samplelist.RData') # "cohort1" "cohort_INDP2" "cohort_INDP3"

load('../Human/04.mndna_model/results/2_FeatureSelection/CRCmodel_fs.RData')
trnval_crc_samplelist = cohort1[which(((cohort1$label=='HD')|(cohort1$label=='CRC'))&(cohort1$group.label=='train/val')), 'Sample']

trnval_crc_df = as.data.frame(t(CRCmodel_fs[, trnval_crc_samplelist]))
trnval_crc_df$label = 'HD'; trnval_crc_df[grep('GLRCRC', rownames(trnval_crc_df_t)), 'label'] = 'CRC'

candidate_fea = read.table('../Human/03.features_annotation/result/GenomicAnnotationIn_CRC0103_addtoRegion.log', head=TRUE, sep='\t', row.names=1)
final_fea = read.table('../Human/04.mndna_model/results/3_FeatureReduction/manu_crc_240103.bed', sep='\t', head=FALSE)
final_feas = str_c('chr',final_fea[,1],':',final_fea[,2],'-',final_fea[,3])
final_feas = candidate_fea[final_feas, c('overlappedGeneNumDensity', 'genename')]

trnval_crc_df = trnval_crc_df[, c(rownames(final_feas), 'label')]
group_means <- aggregate(. ~ label, data = trnval_crc_df, FUN = mean)
mean_diff <- group_means[2:ncol(group_means)][group_means$label == "CRC", ] - group_means[2:ncol(group_means)][group_means$label == "HD", ]

mean_diff <- t(mean_diff)
final_feas = merge(final_feas, mean_diff, by='row.names')
final_feas = final_feas[order(final_feas[,"1"]), ,drop=FALSE]
final_feas$gene_label = str_c(str_extract(final_feas$genename, "^[^,]+"), ',...')
final_feas[is.na(final_feas$gene_label), 'gene_label'] = ''
final_feas[is.na(final_feas$overlappedGeneNumDensity), 'overlappedGeneNumDensity'] = ''
final_feas[final_feas$overlappedGeneNumDensity==1, 'gene_label'] = gsub(',...','',final_feas[final_feas$overlappedGeneNumDensity==1, 'gene_label'])
final_feas$description = str_c(final_feas$Row.names,'\n(',final_feas$overlappedGeneNumDensity,'genes: ', final_feas$gene_label,')')
final_feas$description = gsub('\\(genes: \\)','',final_feas$description)

trnval_crc_df_melt = melt(trnval_crc_df, id.vars=c('label'))
trnval_crc_df_melt$label = factor(trnval_crc_df_melt$label, levels=c('HD','CRC'))
trnval_crc_df_melt = merge(trnval_crc_df_melt, final_feas, by.x='variable', by.y='Row.names')
trnval_crc_df_melt$variable = factor(trnval_crc_df_melt$variable, levels = final_feas$Row.names)
trnval_crc_df_melt = trnval_crc_df_melt[order(trnval_crc_df_melt$variable), ]
trnval_crc_df_melt$description = factor(trnval_crc_df_melt$description, levels = unique(trnval_crc_df_melt$description))
p <- ggplot(trnval_crc_df_melt, aes(x=label,y=value,fill=label))+geom_boxplot(outlier.shape=NA)+
      facet_wrap(~description, scale='free', ncol=8)+scale_fill_manual(values=c('#374E55FF','#E64B35FF'))+ylab('Normalized read counts')+
      stat_compare_means(comparisons=list(c('HD','CRC')), label='p.signif',method='wilcox.test')+
      theme_bw()+theme(strip.text = element_text(size = 5), strip.background=element_blank(), axis.text=element_text(size=5), axis.title.y=element_text(size=8), axis.title.x=element_blank(), legend.title=element_blank(), legend.text=element_text(size=8), panel.grid.major=element_blank(), panel.grid.minor=element_blank())
ggsave(str_c('./FigS6/FigS6b.pdf'),p,width=10, height=5)


# pvalue = c()
# symp = c()
# mean_all = c(); median_all = c(); sd_all = c()
# for(i in c(colnames(trnval_crc_df)[grep('chr',colnames(trnval_crc_df))])){
#   tmp = trnval_crc_df[, c(i, 'label')]; colnames(tmp)[1] = 'fea'
#   tmp$label = factor(tmp$label)
#   pvalue = c(pvalue, wilcox.test(fea ~ label, data=tmp)$p.value)
#   symp <- c(symp, symnum(wilcox.test(fea ~ label, data=tmp)$p.value, corr = FALSE, cutpoints = c(0,0.0001, .001,.01,.05, .1, 1), symbols = c("****","***","**","*","ns",".")))
#   avg = aggregate(. ~ label, data = tmp, FUN = mean)
#   mean_all = rbind(mean_all, c(avg[avg$label=='HD','fea'], avg[avg$label=='CRC','fea']))
#   med = aggregate(. ~ label, data = tmp, FUN = median)
#   median_all = rbind(median_all, c(med[med$label=='HD','fea'], med[med$label=='CRC','fea']))
#   sde = aggregate(. ~ label, data = tmp, FUN = sd)
#   sd_all = rbind(sd_all, c(sde[sde$label=='HD','fea'], sde[sde$label=='CRC','fea']))
# }
# mean_all = as.data.frame(mean_all); colnames(mean_all) = c('HD(mean)','CRC(mean)')
# median_all = as.data.frame(median_all); colnames(median_all) = c('HD(median)','CRC(median)')
# sd_all = as.data.frame(sd_all); colnames(sd_all) = c('HD(sd)','CRC(sd)')

# feature = cbind(c(colnames(trnval_crc_df)[grep('chr',colnames(trnval_crc_df))]), cbind(mean_all, median_all, sd_all))
# feature$p.value = pvalue
# feature$symp = symp
# feature[feature$symp=='.','symp'] = 'ns'
# write.xlsx(feature, './FigS6/TableS6(taMN-DNAfeatures)2.xlsx')