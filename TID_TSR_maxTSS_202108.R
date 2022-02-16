setwd("~/Box Sync/Danko_lab_work/F1_8Tissues/transcription_level_analysis/TID_TSR_maxTSS")
Head="LV"

df=read.table(paste(Head,"_maxTSN_TSS_TID_withSNP_temp3.bed",sep=""))
# 1-9 TSS, 10-15 TID, 17 maxTSN with SNP in TSS(A)

dim(df)
# high low determined by maxTSN allelic read counts
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
# determined high low based on maxTSN allelic read counts
temp <- data.frame(do.call(rbind, strsplit(as.character(df_control$V4), ",")))
df_control$winP = temp[,1]

TSS = unique(df_control[,c(12:20,32)])
dim(TSS)
#temp <- data.frame(do.call(rbind, strsplit(as.character(TSS$V15), ",")))
# TSS$winP = temp[,1]
TSS$log2_high_low = log2((TSS$V18 +1 )/(TSS$V19+1))
TSS$log2_high_low[TSS$winP=="P"] = log2((TSS$V18 +1 )/(TSS$V19+1))[TSS$winP=="P"]

# change the "Head" above
df_BN=df
TSS_BN=TSS
df_LV=df
TSS_LV=TSS


pdf(paste("BN_LV", "_TSSinTID_withSNPsinMaxTSNofTSS_boxplot-2.pdf",sep=""), width=7, height = 7, useDingbats=FALSE)
par(mar=c(6.1, 7.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
par(mgp=c(3,1,0))
par(cex.lab=2.2, cex.axis=2.2)
boxplot(df_BN$log2_high_low, TSS_BN$log2_high_low, df_LV$log2_high_low, TSS_LV$log2_high_low, 
        main=Head, ylab="log2(high allele +1/low allele +1)",
        frame.plot=F,
        outline=FALSE,
        names=c("with SNPs at CA Inr", "control", "LV SNP", "LV control"),
        col=c("purple", "gray"),
        las=2
        )
#abline(h=0, col="red")
abline(h=0, lty=2)
legend("topleft", legend=c("maxTSNs with SNPs", "maxTSNs without SNP" ),
       #title = "SNPs",
       fill = c("purple", "gray"), 
       bty = "n")
dev.off()

length(df_BN$log2_high_low)
length(TSS_BN$log2_high_low)
wilcox.test(df_BN$log2_high_low, TSS_BN$log2_high_low, alternative = "two.sided")

length(df_LV$log2_high_low)
length(TSS_LV$log2_high_low)
wilcox.test(df_LV$log2_high_low, TSS_LV$log2_high_low, alternative = "two.sided")

###########
library(seqLogo)
proportion <- function(x){
        rs <- sum(x);
        return(x / rs);
}
seq_upperCase <- function(seq){
        seq[seq=="a"]<- "A"
        seq[seq=="t"]<- "T"
        seq[seq=="c"]<- "C"
        seq[seq=="g"]<- "G"
        
        a <- NULL
        t <- NULL
        c <- NULL
        g <- NULL
        for (i in 1:NCOL(seq)){
                a <- c(a, sum(seq[,i]=="A"))
                t <- c(t, sum(seq[,i]=="T"))
                c <- c(c, sum(seq[,i]=="C"))
                g <- c(g, sum(seq[,i]=="G"))
        }
        return (data.frame(a,c,g,t))
}

SeqLogo <- function(seq, output) {
        #seq=m$HighAlleleSeq
        seq<- data.frame(do.call(rbind, strsplit(as.character(seq), "")))
        df <- seq_upperCase(seq)
        #create position weight matrix
        pwm <- apply(df, 1, proportion)
        p = makePWM((pwm))
        #p <- makePWM(pwm)
        # slotNames(p)
        # p@consensus
        # p@ic
        # p@width
        # p@alphabet
        pdf(output)
        seqLogo(p)
        dev.off()
        return (pwm)
}

Head="BN"
df=read.table(file = paste(Head, "_maxTSN.TATA_TSS_TID.bed",sep=""), header = F, comment.char = "\"")
#View(df)
colnames(df)[7:9]=c("mat_allele_count","pat_allele_count","identical_reads_count" )
colnames(df)[11:12]=c("B6_seq",  "CAST_seq")
colnames(df)[4] = "hmm_state"
library("TmCalculator")
df$maxTSNlevel = df$mat_allele_count + df$pat_allele_count + df$identical_reads_count
df$B6_AT = unlist(lapply(df$B6_seq, function(x){ 100 - GC(as.character(x))}))
df$CAST_AT = unlist(lapply(df$CAST_seq, function(x){ 100 - GC(as.character(x))}))

df$B6_T = unlist(lapply(df$B6_seq, function(x){ T(as.character(x))}))
df$B6_A = unlist(lapply(df$B6_seq, function(x){ A(as.character(x))}))
df$CAST_A = unlist(lapply(df$CAST_seq, function(x){ A(as.character(x))}))
df$CAST_T = unlist(lapply(df$CAST_seq, function(x){ T(as.character(x))}))
cut=1/8*100
df=df[df$CAST_A > cut | df$CAST_T >cut | df$B6_A > cut | df$B6_T > cut,]

plot(log10(df$maxTSNlevel), log(df$B6_AT))
plot((df$maxTSNlevel)[df$maxTSNlevel>500], (df$B6_AT)[df$maxTSNlevel>500])

SeqLogo(df$B6_seq[df$maxTSNlevel<200], "200-.pdf")
SeqLogo(df$B6_seq[df$maxTSNlevel>100], "100+.pdf")
SeqLogo(df$B6_seq[df$maxTSNlevel>500], "500+.pdf")
SeqLogo(df$B6_seq[df$maxTSNlevel>1500], "1500+.pdf")
SeqLogo(df$B6_seq[df$maxTSNlevel>1000], "1000+.pdf")
SeqLogo(df$B6_seq[df$maxTSNlevel>2000], "2000+.pdf")
sum(df$maxTSNlevel>500)

control = df[(df$B6_AT - df$CAST_AT ==0) ,]
sub_df=df[(df$B6_AT - df$CAST_AT !=0) & (df$mat_allele_count + df$pat_allele_count >0) ,]

#View(sub_df)

# high low determined by AT content of maxTSN upstream
sub_df$log2_high_low = log2((sub_df$mat_allele_count+1 )/(sub_df$pat_allele_count+1))
sub_df$log2_high_low [(sub_df$B6_AT - sub_df$CAST_AT)<0 ] = log2((sub_df$pat_allele_count+1)/ (sub_df$mat_allele_count+1 ))[(sub_df$B6_AT - sub_df$CAST_AT)<0 ]

# TSS high low determined by AT content of maxTSN upstream
sub_df$log2_high_low = log2((sub_df$V19+1 )/(sub_df$V20+1))
sub_df$log2_high_low [(sub_df$B6_AT - sub_df$CAST_AT)<0 ] = log2((sub_df$V20+1)/ (sub_df$V19t+1 ))[(sub_df$B6_AT - sub_df$CAST_AT)<0 ]

# high low determined by maxTSN allelic read counts
temp <- data.frame(do.call(rbind, strsplit(as.character(sub_df$hmm_state), ",")))
sub_df$maxTSN_winP = temp[,1]
sub_df$log2_high_low = log2((sub_df$mat_allele_count+1 )/(sub_df$pat_allele_count+1))
sub_df$log2_high_low [sub_df$maxTSN_winP=="P" ] = log2((sub_df$pat_allele_count+1)/ (sub_df$mat_allele_count+1 ))[sub_df$maxTSN_winP=="P" ]

# TSS high low determined by maxTSN allelic read counts
temp <- data.frame(do.call(rbind, strsplit(as.character(sub_df$hmm_state), ",")))
sub_df$maxTSN_winP = temp[,1]
sub_df$log2_high_low = log2((sub_df$V19+1 )/(sub_df$V20+1))
sub_df$log2_high_low [sub_df$maxTSN_winP=="P" ] =  log2((sub_df$V20+1)/ (sub_df$V19t+1 ))[sub_df$maxTSN_winP=="P" ]



#hist(sub_df$log2_high_low, breaks = seq(-10.25,10.25,0.5) )
 
boxplot(sub_df$log2_high_low, 
        main=Head, ylab="log2(high allele +1/low allele +1)",
        frame.plot=F,
        outline=FALSE,
        #names=c("with SNPs at CA Inr", "control"),
        las=2
)
abline(h=0, col="red")

write.table(sub_df, file=paste(Head,"_maxTSN.TATA_TSS_TID_diffAT.bed", sep = ""),
            quote=F, sep="\t", row.names = F,
            col.names = F)

#TSS within the same TID of the TSS_maxTSN_withDiffTATA
Head="BN"
TSS = read.table(paste(Head,"_maxTSN.TATA_TSS_TID_diffAT_TSSC-A.bed",sep=""))
colnames(TSS)[18:19]=c("B6_seq",  "CAST_seq")
temp <- data.frame(do.call(rbind, strsplit(as.character(TSS$V17), ",")))
TSS$MaxTSNwinP = temp[,1]

library("TmCalculator")
A <- function(a){
        a=s2c(as.character(a))
        return (sum(a=="A"|a=="a")/length(a)*100)
}
T <- function(a){
        a=s2c(as.character(a))
        return (sum(a=="T"|a=="t")/length(a)*100)
}

TSS$B6_AT = unlist(lapply(TSS$B6_seq, function(x){ 100 - GC(as.character(x))}))
TSS$CAST_AT = unlist(lapply(TSS$CAST_seq, function(x){ 100 - GC(as.character(x))}))
TSS$AveAT = (TSS$B6_AT + TSS$CAST_AT)/2

TSS$B6_A = unlist(lapply(TSS$B6_seq, function(x){ A(as.character(x))}))
TSS$B6_T = unlist(lapply(TSS$B6_seq, function(x){ T(as.character(x))}))
TSS$CAST_A = unlist(lapply(TSS$CAST_seq, function(x){ A(as.character(x))}))
TSS$CAST_T = unlist(lapply(TSS$CAST_seq, function(x){ T(as.character(x))}))
dim(TSS)
cut=1/8*100
TSS=TSS[TSS$CAST_A > cut | TSS$CAST_T >cut | TSS$B6_A > cut | TSS$B6_T > cut,]
TSS=TSS[TSS$MaxTSNwinP != "S",]

dim(TSS)
hist(TSS$AveAT)

TSS$log2_high_low = log2((TSS$V7 +1 )/(TSS$V8+1))
TSS$log2_high_low[TSS$MaxTSNwinP=="P"] = log2((TSS$V8 +1 )/(TSS$V7+1))[TSS$MaxTSNwinP=="P"] 
bin=20
vioplot(TSS$log2_high_low, TSS$log2_high_low[TSS$AveAT <bin*1],
        TSS$log2_high_low[ TSS$AveAT >= bin*1 &TSS$AveAT<bin*2 ],
        TSS$log2_high_low[ TSS$AveAT >= bin*2 &TSS$AveAT<bin*3 ],
        TSS$log2_high_low[ TSS$AveAT >= bin*3 &TSS$AveAT<bin*4],
        TSS$log2_high_low[ TSS$AveAT >= bin*4],
        main=Head, ylab="log2(high allele +1/low allele +1)",
        frame.plot=F,
        outline=FALSE,
        #names=c("all", "),
        las=2
)
abline(h=0, col="red")
dim(TSS)[1]
sum(TSS$AveAT <bin*1)
sum(TSS$AveAT >= bin*1 &TSS$AveAT<bin*2)
sum(TSS$AveAT >= bin*2 &TSS$AveAT<bin*3)
sum(TSS$AveAT >= bin*3 &TSS$AveAT<bin*4)
sum(TSS$AveAT >= bin*4)

### 
Head="BN"
df=read.table(paste(Head,"_allReads_TSS_maxTSNs_binomtest_-35To-20INTERSECTmotifM00216_maxScore.bed",sep=""))
colnames(df)[16:19]=c("score.b6","score.cast","strand","score.ave")
colnames(df)[7:9]=c("matRead","patRead", "ideRead")
df$scoreB6minusCAST = df$score.b6 - df$score.cast
df$totalRead = df$matRead + df$patRead + df$ideRead
dim(df)
df=df[df$matRead + df$patRead >0,]
dim(df)
plot((df$matRead - df$patRead)/(df$matRead + df$patRead), df$scoreB6minusCAST)
hist( df$scoreB6minusCAST, breaks = seq(-3,3, 0.1))

plot(((df$matRead - df$patRead)/(df$matRead + df$patRead))[abs(df$scoreB6minusCAST)>0.1], df$scoreB6minusCAST[abs(df$scoreB6minusCAST)>0.1])
cor.test(((df$matRead - df$patRead)/(df$matRead + df$patRead))[abs(df$scoreB6minusCAST)>0.1], df$scoreB6minusCAST[abs(df$scoreB6minusCAST)>0.1])


df=read.table(paste(Head,"_allReads_TSS_maxTSNs_binomtest_-35To-20INTERSECTmotifM09433_maxScore.bed",sep=""))
colnames(df)[16:19]=c("score.b6","score.cast","strand","score.ave")
colnames(df)[7:9]=c("matRead","patRead", "ideRead")
df$scoreB6minusCAST = df$score.b6 - df$score.cast
df$totalRead = df$matRead + df$patRead + df$ideRead
dim(df)
df=df[df$matRead + df$patRead >0,]
dim(df)
plot((df$matRead - df$patRead)/(df$matRead + df$patRead), df$scoreB6minusCAST)
plot(((df$matRead - df$patRead)/(df$matRead + df$patRead))[abs(df$scoreB6minusCAST)>0.1], df$scoreB6minusCAST[abs(df$scoreB6minusCAST)>0.1])
cor.test(((df$matRead - df$patRead)/(df$matRead + df$patRead))[abs(df$scoreB6minusCAST)>0.1], df$scoreB6minusCAST[abs(df$scoreB6minusCAST)>0.1])
