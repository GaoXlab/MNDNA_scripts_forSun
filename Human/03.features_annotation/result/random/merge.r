library(stringr)
order_row = c("weak.enhancers","transcription","enhancers","quescient","transcribed.and.enhancer","polycomb.repressed","acetylations","promoters","weak.transcription","bivalent.promoters","exon","HET","znf","TSS","others","DNase")
merge_all = data.frame(order_row)

random_path = './'
for(i in 0:9){
    df=read.table(str_c(random_path, '/all.',as.character(i),'.grouped_output.txt'),sep='\t',head=T, row.names=1)
    df=df[apply(df,1,sum)<101, ]
    print(dim(df))
    tmp = t(df)[intersect(order_row, rownames(t(df))), ]
    cal_perc = c()
    for(row in order_row){
        cal_perc = c(cal_perc, table(tmp[row,]>50)['TRUE']/dim(tmp)[2])
    }
    cal_perc = as.data.frame(t(t(cal_perc)))
    cal_perc$label = order_row
    # write.table(cal_perc, str_c('all.',as.character(i),'/','cal_perc.log'), sep='\t')
    merge_all = cbind(merge_all, cal_perc)
}
merge_all_tmp = merge_all[grep('V1', colnames(merge_all))]
rownames(merge_all_tmp) = merge_all[,'label']
write.table(merge_all_tmp, str_c(random_path, '/all10times.cal_perc.log'), sep='\t')
