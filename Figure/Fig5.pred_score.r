
## Mouse predictive score
library(dplyr)
library(ggplot2)
library(cowplot)
library(pROC)

trnval = read.table('../Mouse/02.mndna_mouse_svcmodel/results/svc_cv_predictivescore.log', sep='\t', head=TRUE)
trnval <- trnval %>% mutate(target = case_when(grepl("A_Ctrl", X) ~ 0, grepl("A_Apc", X) ~ 1))
trnval <- trnval %>% mutate(label = case_when(grepl("A_Ctrl", X) ~ 'Control (HFD)', grepl("A_Apc", X) ~ 'APC (HFD)'))
roc <- roc(trnval$target, trnval$X0, levels=c(0, 1))
threshold <- coords(roc, "best")['threshold']

test = read.table('../Mouse/02.mndna_mouse_svcmodel/results/svc_test_predictivescore.log', sep='\t', head=TRUE)
test$X = gsub('.8m.nodup.q30|.8m.mm10.nodup.q30','',test$X)
test_sampleinfo = read.table('../Mouse/02.mndna_mouse_svcmodel/modelData/sampleinfo.txt', sep='\t', head=TRUE)
test = merge(test, test_sampleinfo, by.x='X', by.y='Sample'); colnames(test)[ncol(test)] = 'testset_label'

### supplementary fig10
col = c('X','X0','label','testset_label')
df_plot = rbind(trnval[, col], test[test$testset_label=='test_20w', col], test[test$testset_label=='10w', col])
df_plot$X0 = round(as.numeric(df_plot$X0), 4)
df_plot$label = factor(df_plot$label, levels=c('Control (HFD)','APC (HFD)'))
p1 <- ggplot(df_plot[df_plot$testset_label=='cross-validation',], aes(x=label, y=X0, fill=label)) + 
        geom_boxplot(outlier.shape=NA)+ geom_jitter(shape=16, position=position_jitter(0.2), size=0.5)+
        scale_fill_manual(values=c('#9ECAE1','#FB6A4A')) + geom_hline(yintercept=as.numeric(threshold), linetype="dashed", color = "red")+
        theme_bw()+theme(legend.position='none', axis.text.x=element_text(angle=45,  hjust = 1), axis.title.x=element_blank())+ylim(0,1)+ylab('Predictive Score')
p2 <- ggplot(df_plot[df_plot$testset_label=='test_20w',], aes(x=label, y=X0, fill=label)) + 
        geom_boxplot()+ geom_jitter(shape=16, position=position_jitter(0.2), size=0.5)+
        scale_fill_manual(values=c('#9ECAE1','#FB6A4A')) + geom_hline(yintercept=as.numeric(threshold), linetype="dashed", color = "red")+
        theme_bw()+theme(legend.position='none', axis.text.x=element_text(angle=45,  hjust = 1), axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title=element_blank())+ylim(0,1)
p3 <- ggplot(df_plot[df_plot$testset_label=='10w',], aes(x=label, y=X0, fill=label)) + 
        geom_boxplot()+ geom_jitter(shape=16, position=position_jitter(0.2), size=0.5)+
        scale_fill_manual(values=c('#9ECAE1','#FB6A4A')) + geom_hline(yintercept=as.numeric(threshold), linetype="dashed", color = "red")+
        theme_bw()+theme(legend.position='none', axis.text.x=element_text(angle=45,  hjust = 1), axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title=element_blank())+ylim(0,1)
g <- plot_grid(p1, p2, p3,ncol=3, rel_widths=c(1.5,1,1))
ggsave('Fig5/FigS10.mouse_predictivescore.pdf', width=3.4, height=3)


### fig5k
library(tidyr)
df_plot = test[test$testset_label=='neutrolization', col]
df_plot$X0 = round(as.numeric(df_plot$X0), 4)
df_plot = df_plot %>% separate(label, into = c("label", "testset_label"), sep = '\\ \\+\\ ')
df_plot$label = factor(df_plot$label, levels=c('Control (HFD)','APC (HFD)'))
p1 <- ggplot(df_plot[df_plot$testset_label=='PBS, neutrolization',], aes(x=label, y=X0, fill=label)) + 
        geom_boxplot(outlier.shape=NA)+ geom_jitter(shape=16, position=position_jitter(0.2), size=0.5)+
        scale_fill_manual(values=c('#9ECAE1','#FB6A4A')) + geom_hline(yintercept=as.numeric(threshold), linetype="dashed", color = "red")+
        theme_bw()+theme(legend.position='none', axis.text.x=element_text(angle=45,  hjust = 1), axis.title.x=element_blank())+ylim(0,1)+ylab('Predictive Score')
p2 <- ggplot(df_plot[df_plot$testset_label=='IgG, neutrolization',], aes(x=label, y=X0, fill=label)) + 
        geom_boxplot()+ geom_jitter(shape=16, position=position_jitter(0.2), size=0.5)+
        scale_fill_manual(values=c('#9ECAE1','#FB6A4A')) + geom_hline(yintercept=as.numeric(threshold), linetype="dashed", color = "red")+
        theme_bw()+theme(legend.position='none', axis.text.x=element_text(angle=45,  hjust = 1), axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title=element_blank())+ylim(0,1)
p3 <- ggplot(df_plot[df_plot$testset_label=='IL-18, neutrolization',], aes(x=label, y=X0, fill=label)) + 
        geom_boxplot()+ geom_jitter(shape=16, position=position_jitter(0.2), size=0.5)+
        scale_fill_manual(values=c('#9ECAE1','#FB6A4A')) + geom_hline(yintercept=as.numeric(threshold), linetype="dashed", color = "red")+
        theme_bw()+theme(legend.position='none', axis.text.x=element_text(angle=45,  hjust = 1), axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title=element_blank())+ylim(0,1)
g <- plot_grid(p1, p2, p3,ncol=3, rel_widths=c(1.8,1,1))
ggsave('Fig5/Fig5k.mouse_predictivescore.pdf', width=2.3, height=3)


### fig5o.
df_plot = test[test$testset_label=='DSS', col]
df_plot$X0 = round(as.numeric(df_plot$X0), 4)
df_plot$label = factor(df_plot$label, levels=c('untreated','DSS-induced'))
p1 <- ggplot(df_plot, aes(x=label, y=X0, fill=label)) + 
        geom_boxplot(outlier.shape=NA)+ geom_jitter(shape=16, position=position_jitter(0.2), size=0.5)+
        scale_fill_manual(values=c('#9ECAE1','#008A3C')) + geom_hline(yintercept=as.numeric(threshold), linetype="dashed", color = "red")+
        theme_bw()+theme(legend.position='none', axis.text.x=element_text(angle=45,  hjust = 1), axis.title.x=element_blank())+ylim(0,1)+ylab('Predictive Score')
ggsave('Fig5/Fig5o.mouse_predictivescore.pdf', width=1.2, height=3)

### fig5o.
df_plot = test[test$testset_label=='Rpl11', col]
df_plot$X0 = round(as.numeric(df_plot$X0), 4)
df_plot$label = factor(df_plot$label, levels=c('Control (Rpl11 +/+)','Rpl11 haploinsufficiency (+/lox)'))
p1 <- ggplot(df_plot, aes(x=label, y=X0, fill=label)) + 
        geom_boxplot(outlier.shape=NA)+ geom_jitter(shape=16, position=position_jitter(0.2), size=0.5)+
        scale_fill_manual(values=c('#9ECAE1','#4E80BA')) + geom_hline(yintercept=as.numeric(threshold), linetype="dashed", color = "red")+
        theme_bw()+theme(legend.position='none', axis.text.x=element_text(angle=45,  hjust = 1), axis.title.x=element_blank())+ylim(0,1)+ylab('Predictive Score')
ggsave('Fig5/FigS12g.mouse_predictivescore.pdf', width=1.2, height=3)

