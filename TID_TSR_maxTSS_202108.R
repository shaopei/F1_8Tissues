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
#boxplot(df$log2_high_low, main=Head, ylab="log2(high+1/low+1)")
#abline(h=0)
#hist(df$log2_high_low, breaks = seq(-5,5,0.5))


df_control = read.table(paste(Head,"_maxTSN_TSS_TID_OutsideTIDwithSNPmaxTSN.bed",sep=""))
# 1-11 maxTSN, 12-20 TSS, 21-29 TID
dim(df_control)
TSS = unique(df_control[,12:20])
dim(TSS)
temp <- data.frame(do.call(rbind, strsplit(as.character(TSS$V15), ",")))
TSS$winP = temp[,1]
TSS$log2_high_low = log2((TSS$V18 +1 )/(TSS$V19+1))
TSS$log2_high_low[TSS$winP=="P"] = log2((TSS$V18 +1 )/(TSS$V19+1))[TSS$winP=="P"]

pdf(paste(Head, "_TSSinTID_withSNPsinMaxTSNofTSS_boxplot-2.pdf",sep=""), width=7, height = 7, useDingbats=FALSE)
par(mar=c(6.1, 7.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
par(mgp=c(3,1,0))
par(cex.lab=2.2, cex.axis=2.2)
boxplot(df$log2_high_low, TSS$log2_high_low, 
        main=Head, ylab="log2(high allele +1/low allele +1)",
        frame.plot=F,
        #outline=FALSE,
        names=c("with SNPs at CA Inr", "control"),
        las=2
        )
abline(h=0, col="red")
dev.off()
