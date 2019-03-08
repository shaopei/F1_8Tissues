setwd("~/Box Sync/Danko_lab_work/F1_8Tissues/")
df=read.csv("proseq_HT_F1_samples_combined.csv")
View(df)
df$cross = "PB6"
df$cross[df$mousePool %in% c("A","F", "G")]="MB6"
df$Sample_ID=as.factor(paste(df$tissue, df$mousePool, sep="_"))
df$tissue_cross=as.factor(paste(df$tissue, df$cross, sep="_"))
new_level=levels(df$tissue_cross)[c(1,2,5,6,11,12,13,14,3,4,15, 16,9,10,7,8 )]
df$tissue_cross <- factor(df$tissue_cross, levels = new_level)


plot(df$tissue.organ, df$mapped.reads, ylab='mapped read counts')
plot(df$tissue.organ, df$map.sep_raw, ylab='mapped read counts / raw_sep read counts')
plot(df$tissue.organ, df$nodup.dup, ylab='nodup_sep read counts / noadapt_sep read counts')


plot(df$tissue_cross, df$mapped.reads/1000000, ylab='mapped read counts (millions)', col=c("pink", "light blue"), las=2)
plot(df$tissue_cross, df$map.sep_raw, ylab='mapped read counts / raw read counts', col=c("pink", "light blue"), las=2)
plot(df$tissue_cross,100*(1- df$nodup.dup), ylab='1 - nodup_sep read counts / noadapt_sep read counts', col=c("pink", "light blue"), las=2)
plot(df$tissue_cross, df$sep_raw/1000000, ylab='raw read counts (millions)', col=c("pink", "light blue"), las=2)
plot(df$tissue_cross, 100*(df$sep_noadapt.sep_raw), ylab='noadapt/raw read counts (%)', col=c("pink", "light blue"), las=2)
plot(df$tissue_cross, 100*(df$mapping.rate..map.sep.), ylab='mapped reads / nodup_sep read counts (%)', col=c("pink", "light blue"), las=2)





df$Sample_ID=paste(df$tissue.organ, df$mousePool, sep = "_")

df_sub = subset(df, df$tissue %in% c("BN", "LV", "LG"))    
df_sub$Sample_ID <- droplevels(as.factor(df_sub$Sample_ID))
summary(df_sub)
par(mar=c(7.1, 6.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
plot(df_sub$Sample_ID, df_sub$mapped.reads/1000000, las=2, ylab="mapped read counts (Millions)")

d1=aggregate(df$mapped.reads, list(df$tissue), sum)
tapply(df$mapped.reads, df$tissue, sum)

d2=aggregate(cbind(map.sep_raw, nodup.dup) ~ tissue, data = df, mean, na.rm = TRUE)
cbind.data.frame(d1,d2)
