source('Rscript/stat.r')

color_defination <- list.load(str_c(workpath,'/color_defination.json'))
immune.combined <- readRDS(str_c(workpath, '/Final_20w_WTAPC_UMAP_30PCA_0.6Resolution.CellAnno_ClassicMarkers.rds'))
Idents(immune.combined) <- immune.combined$celltype
DefaultAssay(immune.combined) <- "integrated"
HSC_subset <- subset(immune.combined, idents='HSC/MPP')
# DefaultAssay(HSC_subset) <- "integrated"
HSC_subset <- ScaleData(HSC_subset, verbose = FALSE)
HSC_subset <- RunPCA(HSC_subset, npcs = 20, verbose = FALSE)
HSC_subset <- RunUMAP(HSC_subset, reduction = "pca", dims = 1:20)
HSC_subset <- FindNeighbors(HSC_subset, reduction = "pca", dims = 1:20)

HSC_subset <- FindClusters(HSC_subset, resolution = 0.6)

#### 1. cellannotation by classic markers
# 1.1 original clusters:
p1 <- DimPlot(HSC_subset, reduction = "umap", label = TRUE, repel = TRUE, label.size = 6, cols=c(paletteer_d("LaCroixColoR::PeachPear", n = length(table(Idents(immune.combined))), type = "continuous"))) #+ NoLegend()
# ggsave('HSC_Anno1_SeuratClusters_withoutAnnotation.pdf', p1, width=6, height=5)

# 1.2 annotation clusters:
DefaultAssay(HSC_subset) <- "RNA"
Idents(HSC_subset) <- factor(HSC_subset$seurat_clusters, levels=c(1,0,3,2,5,10,11,6,7,8,9,4))
markers = c('Kit','Ly6a','Cd34','Flt3','Fcgr3','Hlf','Shisa5','Fos','Mpo','Elane','S100a8','Ms4a2','Itga2b','Apoe','Car1')
p2 <- DotPlot(HSC_subset, features = markers, col.min=-1, cols=c('white','red','red4'), dot.scale = 3) + RotatedAxis() +theme_bw()+
        theme(axis.text.x=element_text(angle = 90,  hjust = 1, vjust = .5)) + ylab('Cell type')
# ggsave('HSCMPP_Anno2_SeuratClusters_withAnnotation.dotplot.pdf', p1, width=6, height=3.8)

celltype_levels = c('HSC', 'MPP2', 'MPP3', 'MPP4')#,'Niche Cells')
celltype_colorlevels = c('#FECFBAFF', '#C196B7FF','#CE7F5BFF', '#5E81ACFF')
color_defination['HSC_celltype'] <- list(celltype_levels)
color_defination['HSC_celltype_color'] <- list(celltype_colorlevels)

HSC_subset$celltype = factor(HSC_subset$celltype, levels=color_defination$HSC_celltype) #color_defination$HSC_celltype
Idents(HSC_subset) = HSC_subset$celltype
pHSC_cellanno <- DimPlot(HSC_subset, reduction = "umap", label = TRUE, repel=TRUE, label.size = 4, pt.size=0.8,cols=color_defination$HSC_celltype_color)+theme_bw() + NoLegend()
# ggsave('HSC_Anno1_SeuratClusters_withAnnotation.pdf', pHSC_cellanno, width=4, height=3)

# 1.3 statistic summary:
sample_labels = c("W_20w_1","W_20w_2","W_20w_3","W_20w_4","A_20w_1","A_20w_2","A_20w_3","A_20w_4")
sample_groups = c('W_20w','A_20w')
celltypes = c('HSC', 'MPP2', 'MPP3', 'MPP4')
colors = c('#FECFBAFF', '#C196B7FF','#CE7F5BFF', '#5E81ACFF')
p2 <- cellnumber_summary(HSC_subset, sample_labels, sample_groups, celltypes, colors, 'HSC')

g <- plot_grid(pHSC_cellanno,p1,p2,ncol=3,rel_widths=c(1.2,2,1.2))
ggsave('../../Figure/Fig6/FigS14.HSCMPP.pdf', g, width=13, height=3.5)

saveRDS(HSC_subset, file = str_c(workpath,'/Final_20w_WTAPC_UMAP_30PCA_0.6Resolution.CellAnno_ClassicMarkers.HSC_subset.rds'))
list.save(color_defination, file = str_c(workpath,'/color_defination.json'))#"va1_result.txt")
print(str_c("\n===== CellAnnotation is finished! =====\n"))
