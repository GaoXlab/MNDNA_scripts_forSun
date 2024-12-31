library(openxlsx)
library(tableone)
library(stringr)
library(ggpubr)
library(cowplot)
source('./plot_function.r')

setwd('/Users/xingyun/Documents/Data Report/DataReport/For manuscript (原液)/For Manu/FinalCheck202412/Figure/')
## sampleinfo:
load('../Human/04.mndna_model/modelData/clinicalinfo.RData') # "cohort1_clinical" "cohort_INDP"
load('../Human/04.mndna_model/modelData/samplelist.RData') # "cohort1" "cohort_INDP2" "cohort_INDP3"

PanC_model_clinical = as.data.frame(rbind(cohort1_clinical, cohort_INDP))
PanC_model_samplelist = as.data.frame(rbind(cohort1, cohort_INDP2, cohort_INDP3))
PanC_model_sampleinfo = merge(PanC_model_samplelist, PanC_model_clinical, by='Individual'); print(dim(PanC_model_sampleinfo))
PanC_model_sampleinfo$group.label = factor(PanC_model_sampleinfo$group.label, levels=c('train/val','test set', 'independent 2','independent 3'))

## PanC model performance
path="../Human/04.mndna_model/results/4_Classification/"
# PanC_trn_pred = load_prediction_results(str_c(path,'panca/xgbm.predict_result_train'), PanC_model_sampleinfo)
PanC_test_pred = load_prediction_results(str_c(path,'panca/xgbm.predict_result'), PanC_model_sampleinfo)
PanC_test_pred = PanC_test_pred[which(PanC_test_pred$group.label=='test set'), ]
PanC_test_pred$label2 = 'PanC'; PanC_test_pred[PanC_test_pred$label=='HD','label2'] = 'HD'

PanC_test_info = merge(PanC_test_pred[,c('Sample','label','label2','final_prob')], cohort1_clinical, by='Sample')

# Fig4a MN-DNA concentration
colnames(PanC_test_info)[grep('concentration', colnames(PanC_test_info))] = 'MNDNA_concentration'
PanC_test_info$MNDNA_concentration = as.numeric(PanC_test_info$MNDNA_concentration)
p1_MNDNAcon <- ggplot(data = PanC_test_info, aes(x = label2, y = MNDNA_concentration, fill=label2)) + ##
    geom_boxplot(outlier.size = 0.2)+scale_fill_manual(values=c("#374E55FF","#E64B35FF"))+
    theme_classic()+theme(legend.position='none',legend.title=element_blank(),axis.title.x = element_blank())+
    stat_compare_means(aes(label = paste0(..method.., ", ", ..p.signif..)), method = "wilcox.test", label.x.npc = 'left', size=2.5)+ylab('MN-DNA concentration (ng/mL)')

# Fig4b RBC,WBC,PLT,HGB
colnames(PanC_test_info)[grep('RBC|WBC|PLT|HGB', colnames(PanC_test_info))] = c('RBC','WBC','PLT','HGB')
p1_score_RBC_cor <- ggscatter(PanC_test_info[PanC_test_info$label!='HD', ], x = 'RBC', y = 'final_prob', size=1, color = 'black',
          add = "reg.line", conf.int = TRUE, cor.coef = TRUE,cor.coef.size = 2.5,add.params = list(color="#374E55FF", type="dashed"), 
          cor.coeff.args = list(label.x.npc = "left", label.y.npc = "bottom", method = "spearman"), xlab = "RBC", ylab = "MN-DNA predictive score") + ylim(0,1)
p1_score_WBC_cor <- ggscatter(PanC_test_info[PanC_test_info$label!='HD', ], x = 'WBC', y = 'final_prob', size=1, color = 'black',
          add = "reg.line", conf.int = TRUE, cor.coef = TRUE,cor.coef.size = 2.5,add.params = list(color="#374E55FF", type="dashed"), 
          cor.coeff.args = list(label.x.npc = "left", label.y.npc = "bottom", method = "spearman"), xlab = "WBC", ylab = "") + ylim(0,1) + theme(axis.title.y = element_blank())
p1_score_PLT_cor <- ggscatter(PanC_test_info[PanC_test_info$label!='HD', ], x = 'PLT', y = 'final_prob', size=1, color = 'black',
          add = "reg.line", conf.int = TRUE, cor.coef = TRUE,cor.coef.size = 2.5,add.params = list(color="#374E55FF", type="dashed"), 
          cor.coeff.args = list(label.x.npc = "left", label.y.npc = "bottom", method = "spearman"), xlab = "PLT", ylab = "") + ylim(0,1) + theme(axis.title.y = element_blank())
p1_score_HGB_cor <- ggscatter(PanC_test_info[PanC_test_info$label!='HD', ], x = 'HGB', y = 'final_prob', size=1, color = 'black',
          add = "reg.line", conf.int = TRUE, cor.coef = TRUE,cor.coef.size = 2.5,add.params = list(color="#374E55FF", type="dashed"), 
          cor.coeff.args = list(label.x.npc = "left", label.y.npc = "bottom", method = "spearman"), xlab = "HGB", ylab = "") + ylim(0,1) + theme(axis.title.y = element_blank())
panel1 = plot_grid(p1_MNDNAcon, p1_score_RBC_cor, p1_score_WBC_cor, p1_score_PLT_cor, p1_score_HGB_cor, ncol=5 , rel_widths=c(1,1.113,1,1,1))


# Fig4c Gender
p1_score_gender <- ggplot(data = PanC_test_info, aes(x = Gender, y = final_prob)) + ##
  geom_boxplot(outlier.colour = NA)+theme_classic()+geom_jitter(shape=16, position=position_jitter(0.2), size=0.5)+
  theme(legend.position='none',axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), strip.background=element_blank())+#axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+#+ylim(0,1) +
  stat_compare_means(aes(label = paste0(..method.., ", ", ..p.signif..)), method = "wilcox.test", label.x.npc = 'left', size=2.5)+facet_grid(.~label2)+ ylab('MN-DNA predictive score') + xlab('Gender')
# Fig4d Smoking, alcohol
p1_smoking <- ggplot(data = PanC_test_info, aes(x = Smoking.state, y = final_prob)) + ##
  geom_boxplot(outlier.colour = NA)+theme_classic()+geom_jitter(shape=16, position=position_jitter(0.2), size=0.5)+
  theme(legend.position='none',axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), strip.background=element_blank())+#axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+#+ylim(0,1) +
  stat_compare_means(aes(label = paste0(..method.., ", ", ..p.signif..)), method = "kruskal.test", label.x.npc = 'left', size=2.5)+facet_grid(.~label2)+ ylab('MN-DNA predictive score') + xlab('History of smoking')
p2_alcohol <- ggplot(data = PanC_test_info, aes(x = Alcohol.state, y = final_prob)) + ##
  geom_boxplot(outlier.colour = NA)+theme_classic()+geom_jitter(shape=16, position=position_jitter(0.2), size=0.5)+
  theme(legend.position='none',axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), strip.background=element_blank())+#axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+#+ylim(0,1) +
  stat_compare_means(aes(label = paste0(..method.., ", ", ..p.signif..)), method = "kruskal.test", label.x.npc = 'left', size=2.5)+facet_grid(.~label2)+ ylab('MN-DNA predictive score') + xlab('History of alcohol consumption')
panel2 = plot_grid(p1_score_gender, p1_smoking, p2_alcohol, ncol=3 )

# Fig4e Age
p1_score_Age_cor <- ggscatter(PanC_test_info, x = 'Age.at.Collection', y = 'final_prob', size=1, color = 'black',
          add = "reg.line", conf.int = TRUE, cor.coef = TRUE,cor.coef.size = 2.5,add.params = list(color="#374E55FF", type="dashed"), 
          cor.coeff.args = list(label.x.npc = "middle", label.y.npc = "bottom", method = "spearman"), xlab = "Age", ylab = "MN-DNA predictive score") + ylim(0,1)

# Fig4f Lesion info
PanC_test_info$Tumor.lesion = as.numeric(PanC_test_info$Tumor.lesion)
PanC_test_info$Average.size = as.numeric(PanC_test_info$Average.size)
p1_score_LesionNum_cor <- ggscatter(PanC_test_info[PanC_test_info$label!='HD', ], x = 'Tumor.lesion', y = 'final_prob', size=1, color = 'black',
          add = "reg.line", conf.int = TRUE, cor.coef = TRUE,cor.coef.size = 2.5,add.params = list(color="#374E55FF", type="dashed"), 
          cor.coeff.args = list(label.x.npc = "left", label.y.npc = "bottom", method = "spearman"), xlab = "Lesion number", ylab = "MN-DNA predictive score") + ylim(0,1)
p1_score_AverageSize_cor <- ggscatter(PanC_test_info[PanC_test_info$label!='HD', ], x = 'Average.size', y = 'final_prob', size=1, color = 'black',
          add = "reg.line", conf.int = TRUE, cor.coef = TRUE,cor.coef.size = 2.5,add.params = list(color="#374E55FF", type="dashed"),  
          cor.coeff.args = list(label.x.npc = "middle", label.y.npc = "bottom", method = "spearman"), xlab = "Average lesion diameter (cm)", ylab = "MN-DNA predictive score") + ylim(0,1)
panel3 = plot_grid(p1_score_Age_cor, p1_score_LesionNum_cor, p1_score_AverageSize_cor, ncol=3 )

Fig4 = plot_grid(panel1, panel2, panel3, ncol=1, rel_heights=c(1,1.3,1))
dir.create('Fig4/')
ggsave('Fig4/Fig4.clinicalinfo_cor_PanCascore.pdf',Fig4, width=9,height=9)