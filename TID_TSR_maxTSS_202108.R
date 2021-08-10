setwd("~/Box Sync/Danko_lab_work/F1_8Tissues/transcription_level_analysis/TID_TSR_maxTSS")

df=read.table("BN_maxTSN_TSS_TID_withSNP_temp3.bed")
dim(df)
temp <- data.frame(do.call(rbind, strsplit(as.character(df$V17), ",")))
df$maxTSN_winP = temp[,1]
df$log2_high_low = log2((df$V7 +1 )/(df$V8+1))
df$log2_high_low[df$maxTSN_winP=="P"] = log2((df$V8 +1 )/(df$V7+1))[df$maxTSN_winP=="P"] 

library(vioplot)
boxplot(df$log2_high_low)
abline(h=0)
hist(df$log2_high_low, breaks = seq(-5,5,0.5))
