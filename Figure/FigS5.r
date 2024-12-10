library(openxlsx)
library(tableone)
library(ggplot2)
library(ggpubr)

## sampleinfo:
load('../Human/04.mndna_model/modelData/clinicalinfo.RData') # "cohort1_clinical" "cohort_INDP"
load('../Human/04.mndna_model/modelData/samplelist.RData') # "cohort1" "cohort_INDP2" "cohort_INDP3"

## Supplementary Fig S5 related Fig3g (trn, val, test)
sampleinfo_3 = merge(cohort1, cohort1_clinical, by='Individual'); print(dim(sampleinfo_3))

## Supplementary Fig S5 related Fig3g (trn, val, test)
colnames(sampleinfo_3)[grep('RBC|WBC|PLT|HGB', colnames(sampleinfo_3))] = c('RBC','WBC','PLT','HGB')
sampleinfo_3$label = factor(sampleinfo_3$label, levels=c('HD','CRC','TC','BC','GC','LC'))
my_comparisons = list(c('HD','CRC'), c('HD','TC'),c('HD','BC'),c('HD','GC'), c('HD','LC'))
p1_RBC<- ggplot(data = sampleinfo_3, aes(x = label, y = RBC, fill=label)) + ##
  geom_boxplot(outlier.size = 0.2)+theme_classic()+theme(legend.position='none',legend.title=element_blank())+xlab('')+ylab('RBC (×10^12/L)') +
  scale_fill_manual(values=c("#374E55FF","#E64B35FF","#008B4599","#63187999","#00828099","#DF8F44FF"))+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test")
p1_WBC<- ggplot(data = sampleinfo_3, aes(x = label, y = WBC, fill=label)) + ##
  geom_boxplot(outlier.size = 0.2)+theme_classic()+theme(legend.position='none',legend.title=element_blank())+xlab('')+ylab('WBC (×10^9/L)') +
  scale_fill_manual(values=c("#374E55FF","#E64B35FF","#008B4599","#63187999","#00828099","#DF8F44FF"))+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test")
p1_PLT<- ggplot(data = sampleinfo_3, aes(x = label, y = PLT, fill=label)) + ##
  geom_boxplot(outlier.size = 0.2)+theme_classic()+theme(legend.position='none',legend.title=element_blank())+xlab('')+ylab('PLT (×10^9/L)') +
  scale_fill_manual(values=c("#374E55FF","#E64B35FF","#008B4599","#63187999","#00828099","#DF8F44FF"))+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test")
p1_HGB<- ggplot(data = sampleinfo_3, aes(x = label, y = HGB, fill=label)) + ##
  geom_boxplot(outlier.size = 0.2)+theme_classic()+theme(legend.position='none',legend.title=element_blank())+xlab('')+ylab('HGB (g/L)') +
  scale_fill_manual(values=c("#374E55FF","#E64B35FF","#008B4599","#63187999","#00828099","#DF8F44FF"))+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test")

g <- plot_grid(p1_RBC, p1_WBC, p1_PLT, p1_HGB, ncol=4)
dir.create('FigS5')
ggsave('FigS5/FigS5.pdf', g, width=10, height=3.3, family="ArialMT")
