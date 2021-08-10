setwd("~/Box Sync/Danko_lab_work/F1_8Tissues/transcription_level_analysis/TID_TSR_maxTSS")
Head="LV"

df=read.table(paste(Head,"_maxTSN_TSS_TID_withSNP_temp3.bed",sep=""))
# 1-9 TSS, 10-15 TID, 17 maxTSN with SNP in TSS(A)

dim(df)
temp <- data.frame(do.call(rbind, strsplit(as.character(df$V17), ",")))
df$maxTSN_winP = temp[,1]
df$log2_high_low = log2((df$V7 +1 )/(df$V8+1))
df$log2_high_low[df$maxTSN_winP=="P"] = log2((df$V8 +1 )/(df$V7+1))[df$maxTSN_winP=="P"] 

library(vioplot)
boxplot(df$log2_high_low, main=Head, ylab="log2(high+1/low+1)")
abline(h=0)
hist(df$log2_high_low, breaks = seq(-5,5,0.5))


df_control = read.table(paste(Head,"_maxTSN_TSS_TID_OutsideTIDwithSNPmaxTSN.bed",sep=""))
# 1-11 maxTSN, 12-20 TSS, 21-29 TID
dim(df_control)
TSS = unique(df_control[,12:20])
dim(TSS)
temp <- data.frame(do.call(rbind, strsplit(as.character(TSS$V15), ",")))
TSS$winP = temp[,1]
TSS$log2_high_low = log2((TSS$V18 +1 )/(TSS$V19+1))
TSS$log2_high_low[TSS$winP=="P"] = log2((TSS$V18 +1 )/(TSS$V19+1))[TSS$winP=="P"]

boxplot(df$log2_high_low, TSS$log2_high_low, 
        main=Head, ylab="log2(high+1/low+1)",
        names=c("withSNPs at CA", "control")
        )
