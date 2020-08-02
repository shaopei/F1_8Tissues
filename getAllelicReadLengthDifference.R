#R --vanilla --slave --args $(pwd) Tissue mat_name pat_name < XXX.R

args=(commandArgs(TRUE))
setwd(args[1])
Tissue=args[2]
mat_name=args[3]
pat_name=args[4]


mat<-read.table(mat_name)
pat<-read.table(pat_name)
#mat<-read.table("HT_matReads_patReads_TSS_maxTSNs_ratio0.5-2_mat.ReadLength.bed")
#pat<-read.table("HT_matReads_patReads_TSS_maxTSNs_ratio0.5-2_pat.ReadLength.bed")

mat_readLength <- NULL
pat_readLength <- NULL
for (i in 1:dim(mat)[1]){
  mat_readLength = c( mat_readLength, as.numeric(strsplit(as.character(mat$V7[i]), ",")[[1]]))
  pat_readLength = c( pat_readLength, as.numeric(strsplit(as.character(pat$V7[i]), ",")[[1]]))
}


readLength=c(mat_readLength, pat_readLength)
pdf(paste(Tissue,"readLength.pdf", sep = "_"))
par(mfrow=c(3,1))
u = max(readLength)+1
hist(readLength, breaks = seq(-0.5,u,1), main=Tissue)
hist(mat_readLength, breaks = seq(-0.5,u,1))
hist(pat_readLength, breaks = seq(-0.5,u,1))
dev.off()

# identify the real length that are the most abundant (maxPause)
# if there is a tie, the shorter read length is reported
for (i in 1:dim(mat)[1]){
  mat$maxPause_RL[i] = as.numeric(names(sort(table(unlist(strsplit(as.character(mat$V7[i]), ","))),decreasing=TRUE)[1]))
  pat$maxPause_RL[i] = as.numeric(names(sort(table(unlist(strsplit(as.character(pat$V7[i]), ","))),decreasing=TRUE)[1]))
}



combine = mat
colnames(combine)[colnames(combine)=="maxPause_RL"] = "mat_maxPause_RL"
combine$pat_maxPause_RL = pat$maxPause_RL[pat$V2 == combine$V2]

# KS test
for (i in 1:dim(mat)[1]){
  m=as.numeric(strsplit(as.character(mat$V7[i]), ",")[[1]])
  p=as.numeric(strsplit(as.character(pat$V7[i]), ",")[[1]])
  combine$p.value[i] = ks.test(m,p) $ p.value
}
combine$p.value.fdr = p.adjust(combine$p.value, method = "fdr")
combine$deltaDistMaxPause = combine$mat_maxPause_RL - combine$pat_maxPause_RL

#View(combine)
pdf(paste(Tissue,"AllelicReadLengthDifference.pdf", sep = "_"))
u=max(abs(combine$deltaDistMaxPause))+1
hist(abs(combine$deltaDistMaxPause)[combine$p.value.fdr > 0.9] # & combine$deltaDistMaxPause != 0], 
     ,breaks = seq(-0.5,u,1),
     freq = F,
     col = "blue",
     xlab="Difference between allelic read length (highest frequency) from same maxTSN",
     main= paste(Tissue, "Read Length difference", sep=" ")
)
hist(abs(combine$deltaDistMaxPause)[combine$p.value.fdr <= 0.1] #  & combine$deltaDistMaxPause != 0], 
     ,breaks = seq(-0.5,u,1),
     freq = F,
     density = 50, col="red",
     add=T
)
legend("topright", 
       legend = c("FDR > 0.9", "FDR <= 0.1"),
       #pch=c(15,15),
       cex=2, 
       lty=c(0,0),
       #bty="n",
       lwd=1.5, 
       density=c(10000,25),
       angle=c(180,45),
       #angle=45,
       fill=c("blue","red")
       , bty = "n"
)
dev.off()
