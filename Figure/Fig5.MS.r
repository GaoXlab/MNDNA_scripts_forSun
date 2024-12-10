library(openxlsx)
library(ggrepel)
library(cowplot)

data <- read.xlsx('../Mouse/03.MS/Protein_annotation_combine.xlsx', sep='\t')
colnames(data) <- gsub('\t',' ', colnames(data))
data_cytokines_anno <- data[unique(c(grep('cytokine', data$"KEGG pathway"), grep('cytokine', data$"Biological Process"), grep('cytokine', data$"Molecular Function"))), ]

## Group
Tumor = colnames(data_cytokines_anno)[grep('^A_Apc', colnames(data_cytokines_anno))]
WT = colnames(data_cytokines_anno)[grep('^A_Ctrl', colnames(data_cytokines_anno))]

## Remove features where more than half of the values are missing. 
result_all <- data_cytokines_anno[, c('Gene name', Tumor, WT)]
rownames(result_all) <- 1:nrow(result_all)
rm_rows = which(rowSums(is.na(result_all[, c(Tumor, WT)])) / ncol(result_all[, c(Tumor, WT)])>0.5)
result_all <- result_all[setdiff(rownames(result_all), rm_rows), ]

## Differentially expressed proteins, remove NA values
fc <- c(); pvalue <- c()
for(row in 1:nrow(result_all)){
	t <- (result_all[row,Tumor]) ; t=t[!is.na(t)]; t=as.numeric(t)
	w <- (result_all[row,WT]) ; w=w[!is.na(w)]; w=as.numeric(w)
	fc <- c(fc, mean(t)/mean(w))
	pvalue <- c(pvalue, wilcox.test(t,w)$'p.value')
}
result_all$foldchange <- fc
result_all$pvalue <- pvalue


pvalue_cutoff <- 0.05
log2fc_cutoff <- 1
result_all <- result_all[, c('Gene name', 'foldchange', 'pvalue')]; colnames(result_all) <- c('GeneName', 'foldchange', 'pvalue')
result_all$lgfc <- log2(result_all$foldchange)
result_all$regulation = ifelse(result_all$pvalue < pvalue_cutoff & abs(result_all$lgfc) >= log2fc_cutoff, ifelse(result_all$lgfc > log2fc_cutoff ,'Up','Down'), 'NoSig')
result_all$label = ifelse(result_all$pvalue < pvalue_cutoff & abs(result_all$lgfc) >= log2fc_cutoff, as.character(result_all$GeneName), ""); dim(result_all)
if(length(unique(result_all$regulation))==3){
	plot_color <- c("black", "grey", "red4")
} else if(unique(result_all$regulation)[2]=="Up"){
	plot_color <- c("grey", "red4")
} else if(unique(result_all$regulation)[2]=="Down"){
	plot_color <- c("grey", "blue4")
}
result_all$label2=''
result_all[result_all$GeneName=='Il18bp','label2'] = 'Il18bp'
plot_color <- c('grey', 'red4')
p <- ggplot(result_all, aes(x = lgfc, y = -log10(pvalue), colour=label2)) + geom_point(alpha=0.4, size=3.5) + scale_color_manual(values=plot_color) +
	geom_vline(xintercept = c(-log2fc_cutoff, log2fc_cutoff), lty=2, col="dimgrey", lwd=0.5) + geom_hline(yintercept = -log10(pvalue_cutoff), lty=2, col="dimgrey", lwd=0.5) +
    labs(x="log2 (fold change)", y="-log10 (p-value)") + theme_bw() + scale_x_continuous(limits = c(-3,3))+
	theme(plot.title = element_text(hjust = 0.5), legend.position="none", legend.title = element_blank(), axis.text = element_text(colour = 'black')) + 
	geom_text_repel(data = result_all, aes(x = lgfc, y = -log10(pvalue), label = label2), max.overlaps = 40, size = 4.5, box.padding = unit(0.5, "lines"), point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE)

ggsave('Fig5/Fig5e.MS.pdf', p, width=3, height=2.4)