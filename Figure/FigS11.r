
library(stringr)
library(ggplot2)
library(ggalt)
library(readxl)
library(stringr)
library(pheatmap)
library(stringr)
library(GenomicRanges)
library(ggvenn)
library(ggrepel)
source('./plot_function.r')


dir.create('./FigS11')

## crc features 
load('../Human/04.mndna_model/modelData/samplelist.RData') # "cohort1" "cohort_INDP2" "cohort_INDP3"
load('../Human/04.mndna_model/results/2_FeatureSelection/CRCmodel_fs.RData')
trnval_crc_samplelist = cohort1[which(((cohort1$label=='HD')|(cohort1$label=='CRC'))&(cohort1$group.label=='train/val')), 'Sample']
trnval_crc_df = as.data.frame(t(CRCmodel_fs[, trnval_crc_samplelist]))
trnval_crc_df$label = 'HD'; trnval_crc_df[grep('GLRCRC', rownames(trnval_crc_df)), 'label'] = 'CRC'
med = aggregate(. ~ label, data = trnval_crc_df, FUN = median)
trnval_crc_fc = as.data.frame(t(med[med$label=='CRC',setdiff(colnames(med), 'label')]/med[med$label=='HD',setdiff(colnames(med), 'label')]))

crc_pos = rownames(trnval_crc_fc[which(trnval_crc_fc[,1]>1), ,drop=FALSE])
crc_neg = rownames(trnval_crc_fc[which(trnval_crc_fc[,1]<1), ,drop=FALSE])

crc_anno = read.table('../Human/03.features_annotation/result/GenomicAnnotationIn_CRC0103_addtoRegion.log', head=TRUE, sep='\t', row.names=1)
crc_pos_genename = unique(unlist(strsplit(crc_anno[crc_pos, 'genename'], ',')))
crc_neg_genename = unique(unlist(strsplit(crc_anno[crc_neg, 'genename'], ',')))


## apc mouse features
mouse_anno = read.table('../Human/03.features_annotation/result/GenomicAnnotationIn_ApcDataMatrix_SignedFindMarkers_addtoRegion.mouse.log', sep='\t', head=TRUE)
mouse_fs = read.table('../Mouse/02.mndna_mouse_svcmodel/results/DataMatrix_SignedFindMarkers_allsigregions.log', sep='\t', head=TRUE)
apc_pos = mouse_anno[which(mouse_anno$region %in% str_c('chr', rownames(mouse_fs[mouse_fs$Significant.padj.avg_lg2FC=='upsig', ]))), 'genename']
apc_neg = mouse_anno[which(mouse_anno$region %in% str_c('chr', rownames(mouse_fs[mouse_fs$Significant.padj.avg_lg2FC=='downsig', ]))), 'genename']

apc_pos_genename = unique(unique(unlist(strsplit(apc_pos, ','))))
apc_neg_genename = unique(unique(unlist(strsplit(apc_neg, ','))))


## FigS11a
x <- list(
  Human = crc_pos_genename,
  Mouse = toupper(apc_pos_genename)
)
p1 <- ggvenn(x, fill_color = c("white","grey"), stroke_size = 0.5, show_percentage = FALSE, set_name_size = 3) 

x <- list(
  Human = crc_neg_genename,
  Mouse = toupper(apc_neg_genename)
)
p2 <- ggvenn(x, fill_color = c("white","grey"), stroke_size = 0.5, show_percentage = FALSE, set_name_size = 3) 
figs11a <- plot_grid(p1,p2,ncol=1)


## FigS11b
library(Hmisc)
filename='2.ApcDataMatrix_SignedFindMarkers.log'
raw = read.table(str_c('./', filename))

apc_pos_con <- capitalize(tolower(intersect(toupper(apc_pos_genename), crc_pos_genename)))
apc_neg_con <- capitalize(tolower(intersect(toupper(apc_neg_genename), crc_neg_genename)))
apc_con_label <- unique(c(apc_pos_con, apc_neg_con))
apc_con_labeled_regions <- rbind(mouse_anno[grep(paste0(apc_pos_con, collapse='|'), mouse_anno$genename), ], 
                                 mouse_anno[grep(paste0(apc_neg_con, collapse='|'), mouse_anno$genename), ])

all <- raw[, c('avg_lg2FC','one_wayAnova.p_adj','Significant.padj.avg_lg2FC')]; rownames(all) = str_c('chr', rownames(all))
all <- merge(all, apc_con_labeled_regions[, c('region','genename')], by.x='row.names', by.y='region', all.x=TRUE)
all$label <- ''; all[is.na(all$genename), 'genename'] = ''
all[grep('Nell1', all$genename), 'label'] = 'Nell1'; all[grep('Nell1', all$genename), 'genename'] = ''
all[grep('Zranb1', all$genename), 'genename'] = 'Zranb1'
all[grep('chr18:76080001-76090000', all$Row.names), 'genename'] = ''
figs11b <- ggplot(all, aes(x=avg_lg2FC, y=-log10(one_wayAnova.p_adj), color=(Significant.padj.avg_lg2FC))) + 
      geom_point(alpha=0.4, size=0.6) + 
      theme_bw() + xlab("Log2 (fold change)") + ylab("-Log10 (adjusted p-value)") +
      theme(plot.title = element_text(size=15,hjust = 0.5), panel.grid=element_blank(), legend.position='none') + 
      scale_colour_manual(values = c('steelblue','gray','brown')) +
      geom_hline(yintercept = -log10(0.01), lty = 'dashed') +
      geom_text_repel(data = all[all$genename!='', ], aes(label = genename), color='black',
                      size = 3, box.padding = unit(1.5, "lines"),
                      point.padding = unit(0., "lines"),
                      segment.color = "black",
                      show.legend = FALSE, max.overlaps = 100,
                      ylim  = c(2, NA))+
      geom_label_repel(data = all[all$label!='', ], aes(label = label), color='red4', 
                      size = 4, box.padding = unit(0.5, "lines"),
                      point.padding = unit(0.5, "lines"),
                      segment.color = "black",fontface="bold", 
                      show.legend = FALSE, max.overlaps = 100,
                      ylim  = c(2, NA))

g = plot_grid(figs11a, figs11b, ncol=2)
ggsave('./FigS11/FigS11.pdf',g,width=6,height=3)
