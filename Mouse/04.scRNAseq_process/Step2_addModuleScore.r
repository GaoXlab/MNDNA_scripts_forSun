source('Rscript/stat.r')

color_defination <- list.load(str_c(workpath,'/color_defination.json'))
immune.combined <- readRDS(str_c(workpath, '/Final_20w_WTAPC_UMAP_30PCA_0.6Resolution.CellAnno_ClassicMarkers.rds'))

# ### 1. cellcycle in annotation clusters
# color_defination['Phase'] <- list(c("G1","G2M","S"))
# color_defination['Phase_color'] <- list(c("#0E5F98FF","#4CBBAAFF","#C6C7A0FF"))
# s.genes <- cc.genes$s.genes
# g2m.genes <- cc.genes$g2m.genes
# s.genes <- c(capitalize(tolower(cc.genes$s.genes)), 'Cenpu') # Cenpu: Mlf1ip(human)
# g2m.genes <- c(capitalize(tolower(cc.genes$g2m.genes)), 'Pimreg','Jpt1') # Pimreg: Fam64a(Human); Jpt1: Hn1(Human)
# DefaultAssay(immune.combined) <- "RNA"
# immune.combined <- CellCycleScoring(immune.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)

### 2. il18相关(2019 paper, Rex et al.)
library(readxl)
DefaultAssay(immune.combined) <- "RNA"
il18signal_gene <- readxl::read_xlsx('Downloads/12079_2019_544_MOESM6_ESM.xlsx') ## Rex et al. (2019)
il18signal_gene <- gsub('\r\n', '', setdiff(as.data.frame(il18signal_gene)[,1],'Gene name'))
il18signal_gene <- capitalize(tolower(il18signal_gene))
gene <- list()
gene$gene <- il18signal_gene#il18_mediated_signaling_pathway_mm10
immune.combined <- AddModuleScore(object=immune.combined, features=gene, ctrl=100, name='CD_features')
colnames(immune.combined@meta.data)[length(colnames(immune.combined@meta.data))] <- 'il18signal'
immune.combined[['il18signal_stage']] <- ifelse(immune.combined@meta.data[,'il18signal'] > mean(immune.combined@meta.data[,'il18signal']), 'High', 'Low')


pathwayi='il18signal'
DefaultAssay(immune.combined) <- "RNA"
mydata <- FetchData(immune.combined, vars=c('UMAP_1','UMAP_2','sample_group','celltype','il18signal_stage',pathwayi))
colnames(mydata)[dim(mydata)[2]] <- 'pathway_label'
mydata$sample_group <- factor(mydata$sample_group, levels=c('W_20w', 'A_20w'))
p1 <- ggplot(mydata[grep('^HSC/MPP$|^CMP$|^MEP$|^Erythroid$|^CLP$|^B$|^T$|^NK$',mydata$celltype),], aes(x=celltype, y=pathway_label, fill=sample_group)) + 
    geom_violin(trim=TRUE, lwd=0.2, draw_quantiles=0.5) + stat_compare_means(aes(group = sample_group), label = "p.signif", method='wilcox.test') + theme_classic() + 
    scale_fill_manual(values=c('#374E55FF','#E64B35FF')) + 
    theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1), axis.title.x = element_blank(), axis.line.y=element_line(linetype=1,color="black",size=0.3), axis.line.x=element_line(linetype=1,color="black",size=0.3)) + ylab(pathwayi)
p2 <- ggplot(mydata[grep('^GMP$|^Monocyte$|^Macrophage$|^DC$|^Neutrophil$|^Basophil$|^Unknown$',mydata$celltype),], aes(x=celltype, y=pathway_label, fill=sample_group)) + 
    geom_violin(trim=TRUE, lwd=0.2, draw_quantiles=0.5) + stat_compare_means(aes(group = sample_group), label = "p.signif", method='wilcox.test') + theme_classic() + 
    scale_fill_manual(values=c('#374E55FF','#E64B35FF')) + 
    theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1), axis.line.y=element_line(linetype=1,color="black",size=0.3), axis.line.x=element_line(linetype=1,color="black",size=0.3)) + ylab(pathwayi)+xlab('Total BM')
g <- plot_grid(p1,p2,ncol=1)
ggsave('../../Figure/Fig6/Fig6c.AddModuleScore_in_IL18pathway.pdf', g, width=6, height=6)


saveRDS(immune.combined, file = str_c(workpath,'/Final_20w_WTAPC_UMAP_30PCA_0.6Resolution.CellAnno_ClassicMarkers.rds'))
list.save(color_defination, file = str_c(workpath,'/color_defination.json'))#"va1_result.txt")
print(str_c("\n===== CellAnnotation is finished! =====\n"))