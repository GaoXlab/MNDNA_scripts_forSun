#!/usr/bin/env Rscript
source('../Mouse/02.mndna_mouse_svcmodel/scripts/utils.r')
## set the right library paths

dir.create('Fig5/')

## Fig5a. mnDNA-mouse
features = read.table('../Mouse/02.mndna_mouse_svcmodel/results/DataMatrix_SignedFindMarkers_40regions.log',sep='\t',head=TRUE)
samples_type1 <- features[,grep('A_Ctrl', colnames(features))]
samples_type2 <- features[,grep('A_Apc', colnames(features))]
print(str_c("GroupA-sample:",dim(samples_type1)[1],";  GroupB-sample:",dim(samples_type2)[1]))
groups_label = 'Con_APC'
pic_name = str_c('Fig5/Fig5a.pheatmap_Con_APC.pdf')
order_result = pheatmap_f_g2(groups_label, transformat(samples_type1), transformat(samples_type2), pic_name)


## Fig5c
df_apc = read.table('../Human/03.features_annotation/result/Apc.grouped_output.txt',sep='\t',head=T,row.names=1)
df_apc_t = t(df_apc)
a = pheatmap(df_apc[str_c('chr', rownames(order_result)), setdiff(colnames(df_apc),c('BivProm','ReprPC.and.openC'))], cluster_rows=FALSE, cluster_cols=TRUE,
        show_colnames = TRUE, show_rownames = FALSE, gaps_row = 20, color = c(colorRampPalette(colors = c("white","#cb181d","firebrick3"))(250)))
ggsave('Fig5/Fig5c.chromHMM.pdf',a,width=4,height=5)



