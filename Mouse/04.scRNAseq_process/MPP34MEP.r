# .libPathsa('D:/singlecell/R')
library(dplyr)
library(Seurat)
library(stringr)
library(cowplot)
library(viridis)
library(ggpubr)
library(SingleR)
library(scater)
library(metap)
library(paletteer)
library(ggsci)
library(Hmisc)
library(reshape2)
library(plyr)
library(scran)
library(org.Mm.eg.db)
library(openxlsx)
library(ggplot2)
library(ggrepel)
library("rlist")

tf = read.xlsx('Final_20w_WTAPC/TF.xlsx')
tf_cofactors = read.table('Final_20w_WTAPC/Mus_musculus_Cof.txt',head=T,sep='\t')

workpath <- "Final_20w_WTAPC/Final_20w_WTAPC_UMAP_30PCA_0.6Resolution/"
setwd(workpath)

immune.combined <- readRDS(str_c(workpath, '/Final_20w_WTAPC_UMAP_30PCA_0.6Resolution.CellAnno_ClassicMarkers.rds'))
DefaultAssay(immune.combined) <- "RNA"

HSC_subset <- readRDS(str_c(workpath, '/Final_20w_WTAPC_UMAP_30PCA_0.6Resolution.CellAnno_ClassicMarkers.HSC_subset.rds'))
DefaultAssay(HSC_subset) <- "RNA"

    Idents(immune.combined) <- 'celltype'
    HSCMPPMEP <- subset(immune.combined, ident=c('HSC/MPP',"MEP"))
	HSCMPPMEP$celltype2 <- as.character(HSCMPPMEP$celltype)

	MPP3 = intersect(colnames(HSCMPPMEP), rownames(HSC_subset@meta.data[HSC_subset@meta.data$celltype=='MPP3',]))
	MPP4 = intersect(colnames(HSCMPPMEP), rownames(HSC_subset@meta.data[HSC_subset@meta.data$celltype=='MPP4',]))

	HSCMPPMEP$celltype2[MPP3] = "MPP3"
	HSCMPPMEP$celltype2[MPP4] = "MPP4"

	table(HSCMPPMEP$celltype2)
	Idents(HSCMPPMEP) = 'celltype2'

	MPP34MEP = subset(HSCMPPMEP, idents=c('MPP3','MPP4',"MEP"))
  	saveRDS(MPP34MEP,'MPP34MEP.rds')


######## Figure 7
    MPP34MEP = readRDS('MPP34MEP.rds')
	MPP34MEP$sample_group_il18 = str_c(MPP34MEP$sample_group,'_',MPP34MEP$il18signal_stage)
	DefaultAssay(MPP34MEP) <- "RNA"

    Idents(MPP34MEP) <- 'il18signal_stage'

	### High
    IL18hi_MPP34MEP <- subset(MPP34MEP, ident='High')
    DefaultAssay(IL18hi_MPP34MEP) <- "RNA"
    Idents(IL18hi_MPP34MEP) <- 'sample_group'
    IL18hi_MPP34MEP_APCWT_onlyposgenes <- FindMarkers(IL18hi_MPP34MEP, ident.1 = 'A_20w', ident.2 = 'W_20w', verbose = FALSE, only.pos = TRUE)
    IL18hi_MPP34MEP_APCWT_all <- FindMarkers(IL18hi_MPP34MEP, ident.1 = 'A_20w', ident.2 = 'W_20w', verbose = FALSE, only.pos = FALSE, min.pct = 0, logfc.threshold = 0, thresh.test = 0)
    IL18hi_MPP34MEP_APCWT = rownames(IL18hi_MPP34MEP_APCWT_onlyposgenes)
	print(intersect(c(tf$Symbol,tf_cofactors$Symbol), IL18hi_MPP34MEP_APCWT))
	
    ### Low
    IL18lo_MPP34MEP <- subset(MPP34MEP, ident='Low')
    DefaultAssay(IL18lo_MPP34MEP) <- "RNA"
    Idents(IL18lo_MPP34MEP) <- 'sample_group'
    IL18lo_MPP34MEP_APCWT_onlyposgenes <- FindMarkers(IL18lo_MPP34MEP, ident.1 = 'A_20w', ident.2 = 'W_20w', verbose = FALSE, only.pos = TRUE)
    IL18lo_MPP34MEP_APCWT_all <- FindMarkers(IL18lo_MPP34MEP, ident.1 = 'A_20w', ident.2 = 'W_20w', verbose = FALSE, only.pos = FALSE, min.pct = 0, logfc.threshold = 0, thresh.test = 0)
    IL18lo_MPP34MEP_APCWT = rownames(IL18lo_MPP34MEP_APCWT_onlyposgenes)
	print(intersect(c(tf$Symbol,tf_cofactors$Symbol), IL18lo_MPP34MEP_APCWT))




tmp = na.omit(IL18hi_MPP34MEP_APCWT_all[c(tf$Symbol,tf_cofactors$Symbol),])
tmp[tmp$p_val < 0.01 & abs(tmp$avg_log2FC) >= 0.5, ]
tmp$regulation = ifelse(tmp$p_val < 0.01 & abs(tmp$avg_log2FC) >= 0.5, ifelse(tmp$avg_log2FC > 0.5 ,'Up','Down'), 'NoSig')
tmp$label = ifelse(tmp$p_val < 0.01 & (tmp$avg_log2FC) > (0.5), as.character(rownames(tmp)), ""); dim(tmp)
if(length(unique(tmp$regulation))==3){
	plot_color <- c("black", "grey", "red4")
} else if(unique(tmp$regulation)[2]=="Up"){
	plot_color <- c("grey", "red4")
} else if(unique(tmp$regulation)[2]=="Down"){
	plot_color <- c("grey", "blue4")
}
plot_color <- c('#3E60AA',"grey",'#E64B35FF')

p <- ggplot(tmp, aes(x = avg_log2FC, y = -log10(p_val), colour=regulation)) +
	geom_point(alpha=0.5, size=0.4) +
	scale_color_manual(values=c('grey', '#E64B35FF')) +
	# 辅助线
	geom_vline(xintercept = c(-0.5,0.5), lty=2, col="dimgrey", lwd=0.3) +
	geom_hline(yintercept = -log10(0.01), lty=2, col="dimgrey", lwd=0.3) +
	# 坐标轴
	labs(x="log2 (fold change)", y="-log10 (p-value)") + 
	theme_bw() + scale_x_continuous(limits = c(-0.75,0.75))+
    geom_text_repel(data=tmp,  aes(x=avg_log2FC, y= -log10(p_val), label=label), color='grey4', size=3, box.padding = unit(0.3, "lines"),point.padding = unit(0.2, "lines"), segment.color = "grey")+
	# 图例
	theme(plot.title = element_text(hjust = 0.5), legend.position="none", legend.title = element_blank(), axis.text = element_text(colour = 'black'), panel.grid.minor = element_blank(), panel.grid.major = element_blank()) #+ 
# ggsave('IL18hi_MPP34MEP_APCvsWT.TF,cofactors.avg_log2FC0.5.p_val0.01.pdf', p, width=3.2, height=3)
ggsave('../../Figure/Fig7/Fig7a.IL18hi_MPP34MEP.avg_log2FC0.5.p_val0.01.pdf', p, width=3.2, height=2.2)



tmp = na.omit(IL18lo_MPP34MEP_APCWT_all[c(tf$Symbol,tf_cofactors$Symbol),])
# tmp[tmp$p_val < 0.01 & abs(tmp$avg_log2FC) >= 0.3, ]
tmp$regulation = ifelse(tmp$p_val < 0.01 & abs(tmp$avg_log2FC) >= 0.3, ifelse(tmp$avg_log2FC > 0.3 ,'Up','Down'), 'NoSig')
tmp$label = ifelse(tmp$p_val < 0.01 & abs(tmp$avg_log2FC) > (0.3), as.character(rownames(tmp)), ""); dim(tmp)
tmp$label = ifelse(tmp$p_val < 0.00000001, as.character(rownames(tmp)), ""); dim(tmp)
if(length(unique(tmp$regulation))==3){
	plot_color <- c("black", "grey", "red4")
} else if(unique(tmp$regulation)=="NoSig"){
    plot_color <- c("grey")
} else if(unique(tmp$regulation)[2]=="Up"){
	plot_color <- c("grey", "red4")
} else if(unique(tmp$regulation)[2]=="Down"){
	plot_color <- c("grey", "blue4")
} 

p <- ggplot(tmp, aes(x = avg_log2FC, y = -log10(p_val), colour=regulation)) +
	geom_point(alpha=0.5, size=0.4) + #
	scale_color_manual(values='grey') +
	# 辅助线
	geom_vline(xintercept = c(-0.3,0.3), lty=2, col="dimgrey", lwd=0.3) +
	geom_hline(yintercept = -log10(0.01), lty=2, col="dimgrey", lwd=0.3) +
	# 坐标轴
	labs(x="log2 (fold change)", y="-log10 (p-value)") + 
	theme_bw() + scale_x_continuous(limits = c(-0.75,0.75))+
	# 图例
	theme(plot.title = element_text(hjust = 0.5), legend.position="none", legend.title = element_blank(), axis.text = element_text(colour = 'black'), panel.grid.minor = element_blank(), panel.grid.major = element_blank()) #+ 
ggsave('../../Figure/Fig7/Fig7a.IL18lo_MPP34MEP.avg_log2FC0.5.p_val0.01.pdf', p, width=3.2, height=2.2)
