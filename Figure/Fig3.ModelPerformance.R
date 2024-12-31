library(stringr)
library(readxl)
library(openxlsx)
library(cowplot)
source('./plot_function.r')

## sampleinfo:
load('../Human/04.mndna_model/modelData/clinicalinfo.RData') # "cohort1_clinical" "cohort_INDP"
load('../Human/04.mndna_model/modelData/samplelist.RData') # "cohort1" "cohort_INDP2" "cohort_INDP3"

PanC_model_clinical = as.data.frame(rbind(cohort1_clinical, cohort_INDP))
PanC_model_samplelist = as.data.frame(rbind(cohort1, cohort_INDP2, cohort_INDP3))
PanC_model_sampleinfo = merge(PanC_model_samplelist, PanC_model_clinical, by='Individual'); print(dim(PanC_model_sampleinfo))
PanC_model_sampleinfo$group.label = factor(PanC_model_sampleinfo$group.label, levels=c('train/val','test set', 'independent 2','independent 3'))
CRC_model_sampleinfo = PanC_model_sampleinfo[which(PanC_model_sampleinfo$label %in% c('HD','AA','CRC')), ]

dir.create('./Fig3/')
## CRC model performance
path="../Human/04.mndna_model/results/4_Classification/"
CRC_trn_pred = load_prediction_results(str_c(path,'crc/xgbm.predict_result_train'), CRC_model_sampleinfo)
CRC_test_pred = load_prediction_results(str_c(path,'crc/xgbm.predict_result'), CRC_model_sampleinfo)

cutoff = CRC_test_pred[which((CRC_test_pred$Target==0) & (CRC_test_pred$group.label=='test set')), ]
cutoff = as.numeric(cutoff[order(cutoff$final_prob), ][45,'final_prob'])
pdf('./Fig3/Fig3b.pdf',width=4,height=4)
p1 <- plot_auc_95CI(CRC_trn_pred, CRC_test_pred[which(CRC_test_pred$group.label=='train/val'), ], CRC_test_pred[which(CRC_test_pred$group.label=='test set'), ]) 
dev.off()
pdf('./Fig3/Fig3d.pdf',width=4,height=4)
p2 <- plot_auc_95CI_ind(CRC_test_pred[which((CRC_test_pred$group.label=='independent 2')&(CRC_test_pred$label!='AA')), ], 
                        CRC_test_pred[which((CRC_test_pred$group.label=='independent 2')&(CRC_test_pred$label!='CRC')), ])
dev.off()

trn_sen = get_sensitivity_inxspe(CRC_trn_pred, cutoff); trn_sen$classify = 'trn'
val_sen = get_sensitivity_inxspe(CRC_test_pred[which(CRC_test_pred$group.label=='train/val'), ], cutoff); val_sen$classify = 'val'
test_sen = get_sensitivity_inxspe(CRC_test_pred[which(CRC_test_pred$group.label=='test set'), ], cutoff); test_sen$classify = 'test'
summary_sen = as.data.frame(rbind(trn_sen, val_sen, test_sen)); colnames(summary_sen)[1] = 'stages'
summary_sen$classify = factor(summary_sen$classify, levels=c('trn','val','test'))
g1 <- ggplot(summary_sen[summary_sen$stages!='Total', ], aes(x=stages, y=SEN, fill=classify)) + 
      geom_bar(stat="identity", color="black", position=position_dodge(), alpha=0.8) + geom_errorbar(aes(ymin=SEN.low, ymax=SEN.up), width=.2, position=position_dodge(.9)) +
      scale_fill_manual(values = c("darkred", "darkblue", "darkgreen")) + theme_bw() + theme(legend.position='bottom')+ ylab('Sensitivity at 96% specificity')
ggsave('./Fig3/Fig3c.pdf', g1, wi=2.5, height=3)


ind_AAsen = get_AAsensitivity_inxspe(CRC_test_pred[which((CRC_test_pred$group.label=='independent 2')&(CRC_test_pred$label!='CRC')), ], cutoff); ind_AAsen$classify = 'ind2_AA'
ind_CRCsen = get_sensitivity_inxspe(CRC_test_pred[which((CRC_test_pred$group.label=='independent 2')&(CRC_test_pred$label!='AA')), ], cutoff); ind_CRCsen$classify = 'ind2_CRC'
colnames(ind_AAsen)[1] = 'type'; colnames(ind_CRCsen)[1] = 'type'; 
summary_ind = as.data.frame(rbind(ind_AAsen, ind_CRCsen))#, ind_HDspe))
summary_ind$type <- factor(summary_ind$type, levels=c('AA_Total','AA_1','AA_vil','AA_HGD','Total','I','II','III'))
g1 <- ggplot(summary_ind, aes(x=type, y=SEN, fill=type)) + 
      geom_bar(stat="identity", color="black", position=position_dodge(), alpha=0.8) +
      geom_errorbar(aes(ymin=SEN.low, ymax=SEN.up), width=.2, position=position_dodge(.9)) +
      scale_fill_manual(values = c("#E18727FF","#FEE8C8","#FDD49E","#FDBB84","red4","#FC9272", "#FB6A4A","#A50F15")) + theme_bw()+theme(legend.position='none') + ylab('True positive rate')

ind_HDspe = get_HDspecificity_inxspe(CRC_test_pred[which((CRC_test_pred$group.label=='independent 2')&(CRC_test_pred$label!='AA')), ], cutoff); ind_HDspe$classify = 'ind2_HD'
colnames(ind_HDspe)[1] = 'type'; colnames(ind_HDspe)[1] = 'type'; 
ind_HDspe$type <- factor(ind_HDspe$type, levels=c('HD_Total','HD_noinfla','HD_infla'))
g2 <- ggplot(ind_HDspe[ind_HDspe$type!='HD_Total', ], aes(x=type, y=SPE, fill=type)) + 
      geom_bar(stat="identity", color="black", position=position_dodge(), alpha=0.8) +
      geom_errorbar(aes(ymin=SPE.low, ymax=SPE.up), width=.2, position=position_dodge(.9)) +
      scale_fill_manual(values = c("#136BB3","#4391C8")) + theme_bw()+theme(legend.position='none') + ylab('True negative rate')
ggsave('./Fig3/Fig3ef.pdf', plot_grid(g1,g2,ncol=2,rel_widths=c(2.4,1)), wi=6, height=3)

## Independent cohort 3
ind3 = get_sensitivity_inxspe(CRC_test_pred[which(CRC_test_pred$group.label=='independent 3'), ], cutoff); ind3$classify = 'ind3'
ind3_HDspe = get_HDspecificity_inxspe(CRC_test_pred[which(CRC_test_pred$group.label=='independent 3'), ], cutoff); ind3_HDspe$classify = 'ind3_HD'


## PanC model performance
path="../Human/04.mndna_model/results/4_Classification/"
PanC_trn_pred = load_prediction_results(str_c(path,'panca/xgbm.predict_result_train'), PanC_model_sampleinfo)
PanC_test_pred = load_prediction_results(str_c(path,'panca/xgbm.predict_result'), PanC_model_sampleinfo)

pdf('./Fig3/Fig3h.pdf',width=4,height=4)
p1 <- plot_auc_95CI(PanC_trn_pred, PanC_test_pred[which(PanC_test_pred$group.label=='train/val'), ], PanC_test_pred[which(PanC_test_pred$group.label=='test set'), ]) 
dev.off()

cutoff = PanC_test_pred[which((PanC_test_pred$Target==0) & (PanC_test_pred$group.label=='test set')), ]
cutoff = as.numeric(cutoff[order(cutoff$final_prob), ][45,'final_prob'])
panc_test_stage = get_sensitivity_inxspe(PanC_test_pred[which(PanC_test_pred$group.label=='test set'), ], cutoff)
panc_test_type = get_sensitivity_inxspe2(PanC_test_pred[which(PanC_test_pred$group.label=='test set'), ], cutoff)
panc_ind_stage = get_sensitivity_inxspe(PanC_test_pred[which((PanC_test_pred$group.label=='independent 2')|(PanC_test_pred$group.label=='independent 3')), ], cutoff)
panc_ind_type = get_sensitivity_inxspe2(PanC_test_pred[which((PanC_test_pred$group.label=='independent 2')|(PanC_test_pred$group.label=='independent 3')), ], cutoff)

cutoff = PanC_test_pred[which((PanC_test_pred$Target==0) & (PanC_test_pred$group.label=='test set')), ]
cutoff = as.numeric(cutoff[order(cutoff$final_prob), ][46,'final_prob'])
panc_test_stage = get_sensitivity_inxspe(PanC_test_pred[which(PanC_test_pred$group.label=='test set'), ], cutoff)
panc_test_type = get_sensitivity_inxspe2(PanC_test_pred[which(PanC_test_pred$group.label=='test set'), ], cutoff)
