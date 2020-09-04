setwd("~/Box Sync/KD_IGV/2020July/")
file_dir="~/Box Sync/KD_IGV/"
SNP.bw <- paste(file_dir, "P.CAST_M.B6_indelsNsnps_CAST.bam.snp.unfiltered_plus.bw", sep="")


df1=read.table("HT_BothAlleleMaxTSNs_ratio0.5-2_map3_mat.pat.IDEreads_map2ref_DistanceTomaxTSN.bed")
df2=read.table("SK_BothAlleleMaxTSNs_ratio0.5-2_map3_mat.pat.IDEreads_map2ref_DistanceTomaxTSN.bed")
df3=read.table("KD_BothAlleleMaxTSNs_ratio0.5-2_map3_mat.pat.IDEreads_map2ref_DistanceTomaxTSN.bed")

df1$Tissue="HT"
df2$Tissue="SK"
df3$Tissue="KD"
df=rbind.data.frame(df1,df2,df3)  



colnames(df)[7:9]=c("mat_map3", "pat_map3" , "ide_map3" )
#View(df)
for (i in 1:dim(df)[1]){
  # pick maxPause
  # if there is a tie, the shorter read length/mapping position is reported
  reads=as.numeric(unlist(c(strsplit(as.character(df$mat_map3[i]), ","), strsplit(as.character(df$pat_map3[i]), ","), strsplit(as.character(df$ide_map3[i]), ","))))
  df$maxPause_map3AllReads[i] = as.numeric(names(sort(table(reads),decreasing=TRUE)[1]))
  if (!is.na(df$mat_map3[i])){
    df$mat_maxPause_map3[i] = as.numeric(names(sort(table(unlist(strsplit(as.character(df$mat_map3[i]), ","))),decreasing=TRUE)[1]))
  }
  if (!is.na(df$pat_map3[i])){
 df$pat_maxPause_map3[i] = as.numeric(names(sort(table(unlist(strsplit(as.character(df$pat_map3[i]), ","))),decreasing=TRUE)[1]))
 df$earlyPause[i]= min(df$mat_maxPause_map3[i],df$pat_maxPause_map3[i])
 df$latePause[i]= max(df$mat_maxPause_map3[i],df$pat_maxPause_map3[i])
   }

}

# with or without indel 
Tissue_list= c("HT", "SK", "KD")
indel<- NULL
for (i in 1:3){
  indel_i=read.table(paste(Tissue_list[i], "_BothAlleleMaxTSNs_ratio0.5-2_RLmatpat_map3refTomaxTSNmatpat_ClosetIndel.bed", sep = ""))
  indel_i$Tissue=Tissue_list[i]
  indel=rbind.data.frame(indel, indel_i)  
}

indel$V5=111
colnames(indel)[14:15]=c("indel", "dist_Indel_maxTSN")
df_indel <- merge(df, indel, sort=F)
df = df_indel


dim(df)
sum(df$maxPause_map3AllReads == df$pat_maxPause_map3, na.rm = T)
sum(df$maxPause_map3AllReads == df$mat_maxPause_map3, na.rm = T)
sum(df$maxPause_map3AllReads == df$mat_maxPause_map3 & df$maxPause_map3AllReads == df$pat_maxPause_map3, na.rm = T)
sum(df$maxPause_map3AllReads != df$mat_maxPause_map3 & df$maxPause_map3AllReads != df$pat_maxPause_map3, na.rm=T)
#View(df[df$maxPause_map3AllReads != df$mat_maxPause_map3 & df$maxPause_map3AllReads != df$pat_maxPause_map3,])

show.window=20; step=1
bed6 <- df[,1:6]
for (i in 1:NROW(bed6)){
  if(bed6[i,6]=="-") {
    bed6[i,3] <- df[i,3] - df$earlyPause[i] + show.window
    bed6[i,2] <- df[i,2] - df$earlyPause[i] - show.window
  } else {
    bed6[i,2] <- df[i,2] + df$earlyPause[i] - show.window
    bed6[i,3] <- df[i,3] + df$earlyPause[i] + show.window
  }
}


#write.table(bed6, file="test.bed", quote = F, sep="\t", row.names = F, col.names = F)
SNPs.show.window <- read_read_mat_SNPs (SNP.bw, bed6[,c(1:6)], step=step, times=1, use.log=FALSE)
dim(SNPs.show.window)
for (i in 1:NROW(bed6)){
  SNP_distance_earlyPause = which(SNPs.show.window[i,]==1)-show.window-1
  
  df$SNP1_distance_earlyPause[i] = SNP_distance_earlyPause[order(abs(SNP_distance_earlyPause), decreasing = F)[1]]
  df$SNP2_distance_earlyPause[i] = SNP_distance_earlyPause[order(abs(SNP_distance_earlyPause), decreasing = F)[2]]
  df$SNP3_distance_earlyPause[i] = SNP_distance_earlyPause[order(abs(SNP_distance_earlyPause), decreasing = F)[3]]
}





# average pause
## average
for (i in 1:dim(df)[1]){
  df$mat_AvePause_map3[i] = mean(as.numeric(unlist(strsplit(as.character(df$mat_map3[i]), ","))))
  df$pat_AvePause_map3[i] = mean(as.numeric(unlist(strsplit(as.character(df$pat_map3[i]), ","))))
}

# KS test result
df$map3.p.value=1
for (i in 1:dim(df)[1]){
  m=as.numeric(unlist(strsplit(as.character(df$mat_map3[i]), ",")))
  p=as.numeric(unlist(strsplit(as.character(df$pat_map3[i]), ",")))
  if (!is.na(m) & !is.na(p) ){
  df$map3.p.value[i] = ks.test(m,p) $ p.value
  }
}

df$map3.p.value.fdr = p.adjust(df$map3.p.value, method = "fdr")

par(mfrow=c(2,2))
step=1
metaplot.SNPsLocation.aroundMaxPause(df[df$map3.p.value.fdr<=0.1,], name=paste("HSK, maxPause map position, step=", step, sep=""), col="red", step=step)
abline(v=0, col="gray")

#metaplot.SNPsLocation.aroundMaxPause(df[df$map3.p.value.fdr<=0.1,], name="HT, maxPause map position", col="black", step=5, add=TRUE)

metaplot.SNPsLocation.aroundMaxPause(df[df$map3.p.value.fdr>0.9,], col="blue", step=step, add=TRUE)
legend("topright", legend=c(paste("fdr <= 0.1, n= ", sum(df$map3.p.value.fdr<=0.1), sep=""), 
                            paste("fdr >  0.9, n= ", sum(df$map3.p.value.fdr>0.9), sep="")),
       col=c("red", "blue"), bty = "n", lty=1)

sum(df$map3.p.value.fdr<=0.1)
sum(df$map3.p.value.fdr>0.9)

metaplot.SNPsLocation.aroundMaxPause <-function(df, name="", use.sum=FALSE, col="red", show.window = 49, step=1 ,add=FALSE, pch_u=19){
  bed6 <- df[,1:6]
  # maxPause location +- window
  
  for (i in 1:NROW(bed6)){
    if(bed6[i,6]=="-") {
      bed6[i,3] <- df[i,3] - df$maxPause_map3AllReads[i] + show.window
      bed6[i,2] <- df[i,2] - df$maxPause_map3AllReads[i] - show.window
    } else {
      bed6[i,2] <- df[i,2] + df$maxPause_map3AllReads[i] - show.window
      bed6[i,3] <- df[i,3] + df$maxPause_map3AllReads[i] + show.window
    }
  }
  
  
  #write.table(bed6, file="test.bed", quote = F, sep="\t", row.names = F, col.names = F)
  SNPs.show.window <- read_read_mat_SNPs (SNP.bw, bed6[,c(1:6)], step=step, times=1, use.log=FALSE)
  
  
  if (use.sum){
    a = colSums(SNPs.show.window )
  }else{
    a = colMeans(SNPs.show.window )
  }
  
  
  #pdf(metaplot.pdf, width=6, height = 6)
  #par(mar=c(6.1, 7.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
  #par(mgp=c(3,1,0))
  #par(cex.lab=2.2, cex.axis=2.2)
  
  if (step >1){
    y <- a
    x <- seq(-1*(show.window)+step/2, show.window, step)[1:length(y)]
  }else{
    y <- a
    x <- seq(-1*(show.window), show.window, step)[1:length(y)]
  }
  if(add){
    lines(x, y, col=col, type="o", pch=pch_u)
  }else{
    plot(x, y, col=col, 
         main = name,
         ylab= "SNPs counts mean",
         type = "o", 
         pch=pch_u,
         xlab="",
         las=1
         #ylim = c(0,0.5)
         
    )
  }
  
  #dev.off()
  return (list(x, y))
}


metaplot.SNPsLocation.aroundEarlyPause <-function(df, name="", use.sum=FALSE, col="red", show.window = 49, step=1 ,add=FALSE, pch_u=19){
  bed6 <- df[,1:6]
  # maxPause location +- window
  
  for (i in 1:NROW(bed6)){
    if(bed6[i,6]=="-") {
      bed6[i,3] <- df[i,3] - df$earlyPause[i] + show.window
      bed6[i,2] <- df[i,2] - df$earlyPause[i] - show.window
    } else {
      bed6[i,2] <- df[i,2] + df$earlyPause[i] - show.window
      bed6[i,3] <- df[i,3] + df$earlyPause[i] + show.window
    }
  }
  
  
  #write.table(bed6, file="test.bed", quote = F, sep="\t", row.names = F, col.names = F)
  SNPs.show.window <- read_read_mat_SNPs (SNP.bw, bed6[,c(1:6)], step=step, times=1, use.log=FALSE)
  dim(SNPs.show.window)
  for (i in 1:NROW(bed6)){
    SNP_distance_earlyPause = which(SNPs.show.window[i,]==1)-show.window-1
    
    df$SNP1_distance_earlyPause[i] = SNP_distance_earlyPause[order(abs(SNP_distance_earlyPause), decreasing = F)[1]]
  df$SNP2_distance_earlyPause[i] = SNP_distance_earlyPause[order(abs(SNP_distance_earlyPause), decreasing = F)[2]]
  df$SNP3_distance_earlyPause[i] = SNP_distance_earlyPause[order(abs(SNP_distance_earlyPause), decreasing = F)[3]]
  }
  
  
  if (use.sum){
    a = colSums(SNPs.show.window )
  }else{
    a = colMeans(SNPs.show.window )
  }
  
  
  #pdf(metaplot.pdf, width=6, height = 6)
  #par(mar=c(6.1, 7.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
  #par(mgp=c(3,1,0))
  #par(cex.lab=2.2, cex.axis=2.2)
  
  if (step >1){
    y <- a
    x <- seq(-1*(show.window)+step/2, show.window, step)[1:length(y)]
  }else{
    y <- a
    x <- seq(-1*(show.window), show.window, step)[1:length(y)]
  }
  if(add){
    lines(x, y, col=col, type="o", pch=pch_u)
  }else{
    plot(x, y, col=col, 
         main = name,
         ylab= "SNPs counts mean",
         type = "o", 
         xlab="",
         pch=pch_u,
         las=1
         #ylim = c(0,0.5)
         
    )
  }
  
  #dev.off()
  return (list(x, y))
}


metaplot.SNPsLocation.aroundLatePause <-function(df, name="", use.sum=FALSE, col="red", show.window = 49, step=1 ,add=FALSE,pch_u=19){
  bed6 <- df[,1:6]
  # maxPause location +- window
  
  for (i in 1:NROW(bed6)){
    if(bed6[i,6]=="-") {
      bed6[i,3] <- df[i,3] - df$latePause[i] + show.window
      bed6[i,2] <- df[i,2] - df$latePause[i] - show.window
    } else {
      bed6[i,2] <- df[i,2] + df$latePause[i] - show.window
      bed6[i,3] <- df[i,3] + df$latePause[i] + show.window
    }
  }
  
  
  #write.table(bed6, file="test.bed", quote = F, sep="\t", row.names = F, col.names = F)
  SNPs.show.window <- read_read_mat_SNPs (SNP.bw, bed6[,c(1:6)], step=step, times=1, use.log=FALSE)
  
  
  if (use.sum){
    a = colSums(SNPs.show.window )
  }else{
    a = colMeans(SNPs.show.window )
  }
  
  
  #pdf(metaplot.pdf, width=6, height = 6)
  #par(mar=c(6.1, 7.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
  #par(mgp=c(3,1,0))
  #par(cex.lab=2.2, cex.axis=2.2)
  
  if (step >1){
    y <- a
    x <- seq(-1*(show.window)+step/2, show.window, step)[1:length(y)]
  }else{
    y <- a
    x <- seq(-1*(show.window), show.window, step)[1:length(y)]
  }
  if(add){
    lines(x, y, col=col, type="o", pch=pch_u)
  }else{
    plot(x, y, col=col, 
         main = name,
         ylab= "SNPs counts mean",
         type = "o", 
         xlab="",
         pch=pch_u,
         las=1
         #ylim = c(0,0.5)
         
    )
  }
  
  #dev.off()
  return (list(x, y))
}


### early pause with at least bp.apart bp difference
tempFunc <-function(bp.apart=5, upto=10){
step=1


metaplot.SNPsLocation.aroundEarlyPause(df[df$map3.p.value.fdr<=0.1 & abs(df$mat_maxPause_map3 - df$pat_maxPause_map3) >= bp.apart & abs(df$mat_maxPause_map3 - df$pat_maxPause_map3) <= upto,], 
                                       name=paste("HT, step=", step,", early and late at least ", bp.apart," bp apart, upto " ,upto, sep=""), col="red", step=step)
 abline(v=0, col="gray")
# abline(v=-5, col="gray")
# abline(v=-10, col="gray")

metaplot.SNPsLocation.aroundEarlyPause(df[df$map3.p.value.fdr>0.9,], col="blue", step=step, add=TRUE)
# legend("topright", legend=c(paste("fdr <= 0.1, n= ", sum(df$map3.p.value.fdr<=0.1 & abs(df$mat_maxPause_map3 - df$pat_maxPause_map3) >0), sep=""), 
#                             paste("fdr >  0.9, n= ", sum(df$map3.p.value.fdr>0.9), sep="")),
#        col=c("red", "blue"), bty = "n", lty=1, pch=19)
# 

metaplot.SNPsLocation.aroundLatePause(df[df$map3.p.value.fdr<=0.1 & abs(df$mat_maxPause_map3 - df$pat_maxPause_map3) >= bp.apart & abs(df$mat_maxPause_map3 - df$pat_maxPause_map3) <= upto,], 
                                      col="orange", step=step, add=TRUE, pch_u = 12)
metaplot.SNPsLocation.aroundLatePause(df[df$map3.p.value.fdr>0.9,], col="green", step=step, add=TRUE, pch_u = 12)

legend("topright", legend=c(paste("Early Pause, fdr <= 0.1, n= ", sum(df$map3.p.value.fdr<=0.1 & abs(df$mat_maxPause_map3 - df$pat_maxPause_map3) >= bp.apart & abs(df$mat_maxPause_map3 - df$pat_maxPause_map3) <= upto), sep=""), 
                            paste("Early Pause, fdr >  0.9, n= ", sum(df$map3.p.value.fdr>0.9), sep=""),
                            paste("Late Pause, fdr <= 0.1, n= ", sum(df$map3.p.value.fdr<=0.1 & abs(df$mat_maxPause_map3 - df$pat_maxPause_map3) >= bp.apart & abs(df$mat_maxPause_map3 - df$pat_maxPause_map3) <= upto), sep=""),
                            paste("Late Pause, fdr >  0.9, n= ", sum(df$map3.p.value.fdr>0.9), sep="")),
       #col=c("red", "blue", "orange", "green"), 
       bty = "n", lty=1, pch=c(19,19,12,12))
}
tempFunc(1,100)

par(mfrow=c(2,2))
for (i in 1:14){
  tempFunc( i ,100)
  }

tempFunc_Ave <-function(bp.apart=5, upto=10, step=3){
  
  
  
  metaplot.SNPsLocation.aroundEarlyPause(df[df$map3.p.value.fdr<=0.1 & abs(df$mat_AvePause_map3 - df$pat_AvePause_map3) >= bp.apart & abs(df$mat_AvePause_map3 - df$pat_AvePause_map3) <= upto,], 
                                         name=paste("HT, step=", step,", early and late at least ", bp.apart," bp apart, upto " ,upto, sep=""), col="red", step=step)
  abline(v=0, col="gray")
  # abline(v=-5, col="gray")
  # abline(v=-10, col="gray")
  
  metaplot.SNPsLocation.aroundEarlyPause(df[df$map3.p.value.fdr>0.9,], col="blue", step=step, add=TRUE)
  # legend("topright", legend=c(paste("fdr <= 0.1, n= ", sum(df$map3.p.value.fdr<=0.1 & abs(df$mat_AvePause_map3 - df$pat_AvePause_map3) >0), sep=""), 
  #                             paste("fdr >  0.9, n= ", sum(df$map3.p.value.fdr>0.9), sep="")),
  #        col=c("red", "blue"), bty = "n", lty=1, pch=19)
  # 
  
  metaplot.SNPsLocation.aroundLatePause(df[df$map3.p.value.fdr<=0.1 & abs(df$mat_AvePause_map3 - df$pat_AvePause_map3) >= bp.apart & abs(df$mat_AvePause_map3 - df$pat_AvePause_map3) <= upto,], 
                                        col="orange", step=step, add=TRUE, pch_u = 12)
  metaplot.SNPsLocation.aroundLatePause(df[df$map3.p.value.fdr>0.9,], col="green", step=step, add=TRUE, pch_u = 12)
  
  legend("topright", legend=c(paste("Early Pause, fdr <= 0.1, n= ", sum(df$map3.p.value.fdr<=0.1 & abs(df$mat_AvePause_map3 - df$pat_AvePause_map3) >= bp.apart & abs(df$mat_AvePause_map3 - df$pat_AvePause_map3) <= upto), sep=""), 
                              paste("Early Pause, fdr >  0.9, n= ", sum(df$map3.p.value.fdr>0.9), sep=""),
                              paste("Late Pause, fdr <= 0.1, n= ", sum(df$map3.p.value.fdr<=0.1 & abs(df$mat_AvePause_map3 - df$pat_AvePause_map3) >= bp.apart & abs(df$mat_AvePause_map3 - df$pat_AvePause_map3) <= upto), sep=""),
                              paste("Late Pause, fdr >  0.9, n= ", sum(df$map3.p.value.fdr>0.9), sep="")),
         col=c("red", "blue", "orange", "green"), bty = "n", lty=1, pch=c(19,19,12,12))
}
tempFunc(5,100)
par(mfrow=c(3,3))
for (i in 1:9){
tempFunc_Ave(i,i+1, step=1)}


abline(v=-10, col="gray")

df_5bp =df[df$map3.p.value.fdr<=0.1 & abs(df$mat_maxPause_map3 - df$pat_maxPause_map3) >= 5,]

Tissue="HT"
indel=read.table(paste(Tissue, "_BothAlleleMaxTSNs_ratio0.5-2_RLmatpat_map3refTomaxTSNmatpat_ClosetIndel.bed", sep = ""))

library("sqldf")
df_5bp_indel <- sqldf ("select df_5bp.V1,df_5bp.V2,df_5bp.V3,df_5bp.V4,df_5bp.V5,df_5bp.V6,earlyPause,latePause,indel.V14,indel.V15 from df_5bp left join indel on df_5bp.V2=indel.V2" )

df_5bp_indel$indel= df_5bp_indel$V15<=df_5bp_indel$latePause
sum(df_5bp_indel$indel)

for (i in 1:dim(df_5bp_indel)[1]){
df_5bp_indel$mat_Idel_length[i]= length(unlist(strsplit(unlist(strsplit(as.character(df_5bp_indel$V14[i]), ","))[1], "")))
df_5bp_indel$pat_Idel_length[i]= length(unlist(strsplit(unlist(strsplit(as.character(df_5bp_indel$V14[i]), ","))[2], "")))
}
df_5bp_indel$Indel_Len = abs(df_5bp_indel$mat_Idel_length - df_5bp_indel$pat_Idel_length)
df_5bp_indel$pasue_dist=df_5bp_indel$latePause - df_5bp_indel$earlyPause

write.table(df_5bp_indel, "temp.txt", quote = F, row.names = F, sep="\t")

par(mfrow=c(2,1))
sub_df=df[df$map3.p.value.fdr<=0.1,]
hist(abs(sub_df$mat_AvePause_map3 - sub_df$pat_AvePause_map3),
     freq = F, 
     #ylim=c(0,0.4),
     las=1,
     breaks = seq(-0.5,50,1),
     col="black", density = 45, main="map3TomaxTSNs KS FDR <= 0.1")

hist(abs(sub_df$mat_maxPause_map3 - sub_df$pat_maxPause_map3),
     freq=F,
     breaks = seq(-0.5,50,1), col="blue", add=T)

hist(abs(sub_df$mat_AvePause_map3 - sub_df$pat_AvePause_map3),
     freq = F, 
     #ylim=c(0,0.4),
     las=1,
     breaks = seq(-0.5,50,1),
     col="black", density = 45, main="map3TomaxTSNs KS FDR <= 0.1", add=T)

legend("topright", 
       legend = c("AvePause", "MaxPause"),
       title = paste("n = ", dim(sub_df)[1], sep=""),
       #pch=c(15,15),
       #cex=2, 
       lty=c(0,0),
       #bty="n",
       lwd=1.5, 
       density=c(45, 10000),
       angle=c(180,45),
       #angle=45,
       fill=c("black","blue")
       , bty = "n"
)

# map3TomaxTSNs KS FDR <= 0.1 AND SNPs at early pause
if(0){
  sub_df_0.9=df[df$map3.p.value.fdr>0.9,]
  hist(abs(sub_df_0.9$mat_AvePause_map3 - sub_df_0.9$pat_AvePause_map3),
       freq = F, 
       #     ylim=c(0,0.5), 
       las=1,
       breaks = seq(-0.5,50,1),main="maxPause map3TomaxTSNs KS FDR > 0.9")
  
sub_df=df[df$map3.p.value.fdr<=0.1 & df$SNP1_distance_earlyPause==0,]
hist(abs(sub_df$mat_AvePause_map3 - sub_df$pat_AvePause_map3),
          freq = F, ylim=c(0,0.4), las=1,
    breaks = seq(-0.5,50,1),col="blue", main="map3TomaxTSNs KS FDR <= 0.1 AND SNPs at early pause")

hist(abs(sub_df$mat_maxPause_map3 - sub_df$pat_maxPause_map3),
     freq=F,
          breaks = seq(-0.5,50,1), col="red", density = 45, add=T)

legend("topright", 
       legend = c("AvePause", "MaxPause"),
       title = paste("n = ", dim(sub_df)[1], sep=""),
       #pch=c(15,15),
       #cex=2, 
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
}

par(mfrow=c(2,1))
sub_df_0.9=df[df$map3.p.value.fdr>0.9,]
hist(abs(sub_df_0.9$mat_maxPause_map3 - sub_df_0.9$pat_maxPause_map3),
     freq = F, 
#     ylim=c(0,0.5), 
     las=1,
     breaks = seq(-0.5,50,1),main="maxPause map3TomaxTSNs KS FDR > 0.9")

sub_df=df[df$map3.p.value.fdr<=0.1,]
hist(abs(sub_df$mat_maxPause_map3 - sub_df$pat_maxPause_map3),
     freq = F, ylim=c(0,0.4), las=1,
     breaks = seq(-0.5,50,1), col="blue", main="maxPause map3TomaxTSNs KS FDR <= 0.1")

sub_df_SNP=df[df$map3.p.value.fdr<=0.1 & df$SNP1_distance_earlyPause==0 & !is.na(df$SNP1_distance_earlyPause),]
hist(abs(sub_df_SNP$mat_maxPause_map3 - sub_df_SNP$pat_maxPause_map3),
     freq = F, ylim=c(0,0.4), las=1,
     breaks = seq(-0.5,50,1),
     col="red", density = 45, add=T)

legend("topright", 
       legend = c( paste("All FDR<=0.1, n = ", dim(sub_df)[1], sep=""), paste("SNP at early pause, n=", dim(sub_df_SNP)[1])),
       title = ,
       #pch=c(15,15),
       #cex=2, 
       lty=c(0,0),
       #bty="n",
       lwd=1.5, 
       density=c(10000,25),
       angle=c(180,45),
       #angle=45,
       fill=c("blue","red")
       , bty = "n"
)
# exclude difference 0
df$AllelicMaxPauseDist = abs(df$mat_maxPause_map3 - df$pat_maxPause_map3)
sub_df=df[df$map3.p.value.fdr<=0.1 & df$AllelicMaxPauseDist >0,]
sub_df_SNP=df[df$map3.p.value.fdr<=0.1 & df$AllelicMaxPauseDist >0 & df$SNP1_distance_earlyPause==0 & !is.na(df$SNP1_distance_earlyPause),]
par(mfrow=c(2,1))
plot(ecdf(sub_df$AllelicMaxPauseDist), col="blue", 
     xlab="AllelicMaxPauseDist",
     ylab="density",
     las=1,
     main="maxPause map3TomaxTSNs KS FDR <= 0.1, e and l >1bp move")
hist(sub_df$AllelicMaxPauseDist,
     freq = F, ylim=c(0,0.4), las=1,
     breaks = seq(-0.5,50,1), col="blue", main="maxPause map3TomaxTSNs KS FDR <= 0.1",
     add=T
     )

hist(sub_df_SNP$AllelicMaxPauseDist,     freq = F, ylim=c(0,0.4), las=1,
     breaks = seq(-0.5,50,1),
     col="red", density = 45, add=T)

lines(ecdf(sub_df$AllelicMaxPauseDist), col="blue")
lines(ecdf(sub_df_SNP$AllelicMaxPauseDist), col="red")


legend("right", 
       legend = c( paste("All FDR<=0.1, >1bp move, n = ", dim(sub_df)[1], sep=""), 
                   paste("SNP at early pause, n=", dim(sub_df_SNP)[1], sep="")),
       title = ,
       #pch=c(15,15),
       #cex=2, 
       lty=c(0,0),
       #bty="n",
       lwd=1.5, 
       density=c(10000,25),
       angle=c(180,45),
       #angle=45,
       fill=c("blue","red")
       , bty = "n"
)


sub_df_indel=df[df$map3.p.value.fdr<=0.1 & df$AllelicMaxPauseDist >0 & df$dist_Indel_maxTSN <= df$latePause,]
sub_df_no_indel=df[df$map3.p.value.fdr<=0.1 & df$AllelicMaxPauseDist >0 & df$dist_Indel_maxTSN > df$latePause,]
sub_df_SNP_indel=df[df$map3.p.value.fdr<=0.1 & df$AllelicMaxPauseDist >0 & df$SNP1_distance_earlyPause==0 & !is.na(df$SNP1_distance_earlyPause)& df$dist_Indel_maxTSN <= df$latePause,]
dim(sub_df_indel)
dim(sub_df_SNP_indel)

plot(ecdf(sub_df$AllelicMaxPauseDist), col="blue", 
     xlab="AllelicMaxPauseDist",
     ylab="density",
     las=1,
     main="maxPause map3TomaxTSNs KS FDR <= 0.1, e and l >1bp move")

hist(sub_df_no_indel$AllelicMaxPauseDist,
     freq = F, ylim=c(0,0.4), las=1,
     breaks = seq(-0.5,50,1), col="dark orange", main="maxPause map3TomaxTSNs KS FDR <= 0.1",
     add=T
)

hist(sub_df_indel$AllelicMaxPauseDist,     freq = F, ylim=c(0,0.4), las=1,
     breaks = seq(-0.5,50,1),
     col="dark green", density = 45, add=T)

lines(ecdf(sub_df$AllelicMaxPauseDist), col="blue")
lines(ecdf(sub_df_indel$AllelicMaxPauseDist), col="dark green")
lines(ecdf(sub_df_no_indel$AllelicMaxPauseDist), col="dark orange")

legend("right", 
       legend = c( paste("All FDR<=0.1, >1bp move, n = ", dim(sub_df)[1], sep=""), 
                   paste("Without Indel, n=", dim(sub_df_no_indel)[1], sep=""),
                   paste("With Indel, n=", dim(sub_df_indel)[1], sep="")),
       title = ,
       #pch=c(15,15),
       #cex=2, 
       #lty=c(0,0),
       #bty="n",
       lwd=1.5, 
       #density=c(10000,25),
       #angle=c(180,45),
       #angle=45,
       #fill=c("blue","dark organe","dark green")
       col=c("blue","dark orange","dark green"),
       pch=16,
       bty = "n"
)

# KS test 
ks.test(sub_df$AllelicMaxPauseDist ,sub_df_SNP$AllelicMaxPauseDist, alternative = "less")
ks.test(sub_df$AllelicMaxPauseDist ,sub_df_indel$AllelicMaxPauseDist, alternative = "greater")
ks.test(sub_df_no_indel$AllelicMaxPauseDist ,sub_df_indel$AllelicMaxPauseDist, alternative = "greater")


# SNPs at -3 G
# exclude difference 0
df$AllelicMaxPauseDist = abs(df$mat_maxPause_map3 - df$pat_maxPause_map3)
sub_df=df[df$map3.p.value.fdr<=0.1 & df$AllelicMaxPauseDist >0,]
sub_df_SNP_3G=df[df$map3.p.value.fdr<=0.1 & df$AllelicMaxPauseDist >0 &((df$SNP1_distance_earlyPause==-3 & !is.na(df$SNP1_distance_earlyPause))|(df$SNP2_distance_earlyPause==-3 & !is.na(df$SNP2_distance_earlyPause))),]
par(mfrow=c(2,1))
plot(ecdf(sub_df$AllelicMaxPauseDist), col="blue", 
     xlab="AllelicMaxPauseDist",
     ylab="density",
     las=1,
     main="maxPause map3TomaxTSNs KS FDR <= 0.1, e and l >1bp move")
hist(sub_df$AllelicMaxPauseDist,
     freq = F, ylim=c(0,0.4), las=1,
     breaks = seq(-0.5,50,1), col="blue", main="maxPause map3TomaxTSNs KS FDR <= 0.1",
     add=T
)

hist(sub_df_SNP_3G$AllelicMaxPauseDist,     freq = F, ylim=c(0,0.4), las=1,
     breaks = seq(-0.5,50,1),
     col="dark orange", density = 45, add=T)

lines(ecdf(sub_df$AllelicMaxPauseDist), col="blue")
lines(ecdf(sub_df_SNP_3G$AllelicMaxPauseDist), col="dark orange")

legend("right", 
       legend = c( paste("SNP at early pause, n=", dim(sub_df_SNP)[1], sep=""), 
                   paste("SNP at -3 early pause, n=", dim(sub_df_SNP_3G)[1], sep=""),
                   paste("All FDR<=0.1, >1bp move, n = ", dim(sub_df)[1], sep="")),
       title = ,
       #pch=c(15,15),
       #cex=2, 
       #lty=c(0,0),
       #bty="n",
       lwd=1.5, 
       #density=c(10000,25),
       #angle=c(180,45),
       #angle=45,
       #fill=c("blue","dark orange","dark green")
       col=c("red","dark orange","blue"),
       pch=16,
       bty = "n"
)


legend("right", 
       legend = c( paste("All FDR<=0.1, >1bp move, n = ", dim(sub_df)[1], sep=""), 
                   paste("SNP at -3 early pause, n=", dim(sub_df_SNP_3G)[1], sep="")),
       title = ,
       #pch=c(15,15),
       #cex=2, 
       lty=c(0,0),
       #bty="n",
       lwd=1.5, 
       density=c(10000,25),
       angle=c(180,45),
       #angle=45,
       fill=c("blue","red")
       , bty = "n"
)
ks.test(sub_df$AllelicMaxPauseDist ,sub_df_SNP_3G$AllelicMaxPauseDist, alternative = "less")

###



# SNPs at other positions
snp_pause_distance<-function(snp_position=-3){
sub_df=df[df$map3.p.value.fdr<=0.1,]
hist(abs(sub_df$mat_maxPause_map3 - sub_df$pat_maxPause_map3),
     freq = F, 
     ylim=c(0,0.5), 
     las=1,
     breaks = seq(-0.5,50,1), col="blue", main="maxPause map3TomaxTSNs KS FDR <= 0.1")


sub_df_SNP=df[df$map3.p.value.fdr<=0.1 & ((!is.na(df$SNP1_distance_earlyPause) & df$SNP1_distance_earlyPause==snp_position) | (!is.na(df$SNP2_distance_earlyPause) & df$SNP2_distance_earlyPause==snp_position)), ]
hist(abs(sub_df_SNP$mat_maxPause_map3 - sub_df_SNP$pat_maxPause_map3),
     freq = F, 
     #ylim=c(0,0.4), 
     las=1,
     breaks = seq(-0.5,50,1),
     col="red", density = 45, add=T)

legend("topright", 
       legend = c( paste("All, n = ", dim(sub_df)[1], sep=""), paste("SNP at", snp_position  , " ,n=", dim(sub_df_SNP)[1])),
       title = ,
       #pch=c(15,15),
       #cex=2, 
       lty=c(0,0),
       #bty="n",
       lwd=1.5, 
       density=c(10000,25),
       angle=c(180,45),
       #angle=45,
       fill=c("blue","red")
       , bty = "n"
)
}
par(mfrow=c(2,3))
for (i in 0:11){
snp_pause_distance(i)
}


### What kind of SNPs?

df$earlyPause_parent = "M"
df$earlyPause_parent[df$earlyPause==df$pat_maxPause_map3] = "P"
df$earlyPause_parent[df$earlyPause==df$mat_maxPause_map3] = "M"
sum(df$earlyPause_parent=="M")
show.window=0
bed7 <- df[,1:6]
bed7$earlyPause_parent = df$earlyPause_parent
for (i in 1:NROW(bed7)){
  if(bed7[i,6]=="-") {
    bed7[i,3] <- df[i,3] - df$earlyPause[i] + show.window
    bed7[i,2] <- df[i,2] - df$earlyPause[i] - show.window
  } else {
    bed7[i,2] <- df[i,2] + df$earlyPause[i] - show.window
    bed7[i,3] <- df[i,3] + df$earlyPause[i] + show.window
  }
}
# early and late pause with at least 1 bp apart AND map3 KS test fdr<=0.1
bed7 = bed7[df$earlyPause != df$latePause & df$map3.p.value.fdr<=0.1,]
dim(bed7)
dim(df)
write.table(bed7, file="Tissues3_EarlyPause_1bpapart_KSfdr0.1.bed", quote = F, sep="\t", row.names = F, col.names = F)


# early pause, background 
bed0 <- df[,1:6]
bed0$earlyPause_parent = df$earlyPause_parent
show.window=0
for (i in 1:NROW(bed0)){
  if(bed0[i,6]=="-") {
    bed0[i,3] <- df[i,3] - df$earlyPause[i] + show.window
    bed0[i,2] <- df[i,2] - df$earlyPause[i] - show.window
  } else {
    bed0[i,2] <- df[i,2] + df$earlyPause[i] - show.window
    bed0[i,3] <- df[i,3] + df$earlyPause[i] + show.window
  }
}
dim(bed0)
dim(df)
write.table(bed0, file="Tissues3_EarlyPause_BG.bed", quote = F, sep="\t", row.names = F, col.names = F)



bed9 <- df[,1:6]
bed9$earlyPause_parent = df$earlyPause_parent
for (i in 1:NROW(bed9)){
  if(bed9[i,6]=="-") {
    bed9[i,3] <- df[i,3] - df$latePause[i] 
    bed9[i,2] <- df[i,2] - df$latePause[i] 
  } else {
    bed9[i,2] <- df[i,2] + df$latePause[i]  
    bed9[i,3] <- df[i,3] + df$latePause[i]  
  }
}
# early and late pause with at least 1 bp apart
bed9 = bed9[df$earlyPause != df$latePause & df$map3.p.value.fdr<=0.1,]
dim(bed9)
dim(df)
write.table(bed9, file="Tissues3_LatePause_1bpapart_KSfdr0.1.bed", quote = F, sep="\t", row.names = F, col.names = F)





bed8 <- df[,1:6]
for (i in 1:NROW(bed8)){
  if(bed8[i,6]=="-") {
    bed8[i,3] <- df[i,3] - df$maxPause_map3AllReads[i]
    bed8[i,2] <- df[i,2] - df$maxPause_map3AllReads[i]
  } else {
    bed8[i,2] <- df[i,2] + df$maxPause_map3AllReads[i]
    bed8[i,3] <- df[i,3] + df$maxPause_map3AllReads[i]
  }
}

dim(bed8)
bed8 = bed8[!duplicated(bed8$V2),]
dim(bed8)
write.table(bed8, file="combine_maxPause_noduplicate.bed", quote = F, sep="\t", row.names = F, col.names = F)



library("TmCalculator")
library(seqLogo)
#define function that divides the frequency by the row sum i.e. proportions
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

SeqLogo <- function(seq, output, range=NULL) {
  #seq=m$HighAlleleSeq
  seq<- data.frame(do.call(rbind, strsplit(as.character(seq), "")))
  df <- seq_upperCase(seq)
  #create position weight matrix
  pwm <- apply(df, 1, proportion)
  if (!is.null(range)){
    p = makePWM((pwm[,range]))    
  }else{
    p = makePWM((pwm))  
  }

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

seq=read.table("Tissues3_EarlyPause_1bpapart_KSfdr0.1_+-10_Early_LateAlleleSeq.bed")
Tissues3_EarlyPause_1bpapart_KSfdr0.1_early=SeqLogo(seq$V8, "Tissues3_EarlyPause_1bpapart_KSfdr0.1_early.pdf")
Tissues3_EarlyPause_1bpapart_KSfdr0.1_late=SeqLogo(seq$V9, "Tissues3_EarlyPause_1bpapart_KSfdr0.1_late.pdf")

seq=read.table("Tissues3_LatePause_1bpapart_KSfdr0.1_+-10_Early_LateAlleleSeq.bed")
Tissues3_LatePause_1bpapart_KSfdr0.1_early=SeqLogo(seq$V8, "Tissues3_LatePause_1bpapart_KSfdr0.1_early.pdf")
Tissues3_LatePause_1bpapart_KSfdr0.1_late=SeqLogo(seq$V9, "Tissues3_LatePause_1bpapart_KSfdr0.1_late.pdf")

seq=read.table("Tissues3_EarlyPause_BG_+-10_Early_LateAlleleSeq.bed")
Tissues3_EarlyPause_BG_early=SeqLogo(seq$V8, "Tissues3_EarlyPause_BG_early.pdf")
Tissues3_EarlyPause_BG_late=SeqLogo(seq$V9, "Tissues3_EarlyPause_BG_late.pdf")


seq=read.table("combine_maxPause_noduplicate_+-30_mm10_Seq.bed")
combine_maxPause_noduplicate= SeqLogo(seq$V7, "combine_maxPause_noduplicate_+-30_mm10_Se.pdf")

acgt_col=c("dark green", "blue", "orange" , "red")
acgt=c("A","C","G","T")

# combine maxPause 
bin=10
# G content 
b_df=data.frame(t(combine_maxPause_noduplicate))
b_df$group=NA
b_df$group[(30-bin*2+1):(30-bin)]= "1 up up"
b_df$group[(30-bin+1):30]="2 upsteram"
b_df$group[31:(30+bin)]="3 downstream"
boxplot(g~group, data=b_df, main=paste("bin size = ", bin, sep=""))
stripchart(g~group, data=b_df, 
           vertical = TRUE, add = TRUE, method = "jitter", 
           pch=1, 
           col="black",jitter = 0.1, offset = 0, 
           #cex=2
           )
boxplot((g+c)~group, data=b_df)
stripchart((g+c)~group, data=b_df, 
           vertical = TRUE, add = TRUE, method = "jitter", 
           pch=1, 
           col="black",jitter = 0.1, offset = 0, 
           #cex=2
)


C <- function(a){
  return (sum(a=="C"|a=="c")/length(a)*100)
}
G <- function(a){
  return (sum(a=="G"|a=="g")/length(a)*100)
}
par(mfrow=c(1,3))
bin=10
#for (bin in c(10,12)){
seq_df=read.table("combine_maxPause_noduplicate_+-30_mm10_Seq.bed")
range1=(30-bin*2+1):(30-bin)
range2=(30-bin+1):30
range3=31:(30+bin)
for (i in 1:NROW(seq)){
  a=s2c(as.character(seq$V7[i]))
  a1=a[range1]
  a2=a[range2]
  a3=a[range3]
  seq_df$GC_range1[i]= GC(a1)
  seq_df$GC_range2[i]= GC(a2)
  seq_df$GC_range3[i]= GC(a3)
  seq_df$G_range1[i]= G(a1)
  seq_df$G_range2[i]= G(a2)
  seq_df$G_range3[i]= G(a3)
  seq_df$C_range1[i]= C(a1)
  seq_df$C_range2[i]= C(a2)
  seq_df$C_range3[i]= C(a3)
  #sum(a=="C"|a=="c")/length(a)*100
}

library("vioplot")
vioplot(seq_df$G_range1, seq_df$G_range2, seq_df$G_range3, 
        main= paste("bin size = ", bin, sep=""),
        ylab="G content", ylim=c(0,100))
abline(h=40, col="yellow")
vioplot(seq_df$GC_range1, seq_df$GC_range2, seq_df$GC_range3, 
        main= paste("bin size = ", bin, sep=""),
        ylab="GC content", ylim=c(0,100))
abline(h=60, col="green")
vioplot(seq_df$C_range1, seq_df$C_range2, seq_df$C_range3, 
        main= paste("bin size = ", bin, sep=""),
        ylab="C content", ylim=c(0,100))
abline(h=20, col="yellow")

}
# Wilcoxon Rank Sum and Signed Rank Tests
wilcox.test(seq_df$G_range1, seq_df$G_range2, paired = TRUE)
wilcox.test(seq_df$G_range2, seq_df$G_range3, paired = TRUE)
wilcox.test(seq_df$G_range1, seq_df$G_range3, paired = TRUE)

wilcox.test(seq_df$GC_range1, seq_df$GC_range2, paired = TRUE)
wilcox.test(seq_df$GC_range2, seq_df$GC_range3, paired = TRUE)
wilcox.test(seq_df$GC_range1, seq_df$GC_range3, paired = TRUE)

wilcox.test(seq_df$C_range1, seq_df$C_range2, paired = TRUE)
wilcox.test(seq_df$C_range2, seq_df$C_range3, paired = TRUE)
wilcox.test(seq_df$C_range1, seq_df$C_range3, paired = TRUE)







#pdf(paste(organ,"_deltaATCG_single.pdf" ,sep=""), width =7, height = 7)
par(mfcol=c(3,1))
par(mar=c(6.1, 7.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
par(mgp=c(3,1,0))
par(cex.lab=2.2, cex.axis=2.2)
#par(cex=1.5)
pch_u=15
w=10
# background
organ="3 Tissues" ; target=Tissues3_EarlyPause_BG_early - Tissues3_EarlyPause_BG_late
plot(-w:w, target[1,] , col = acgt_col[1], type="o",
     ylim=c(-0.15,0.15), 
     pch=pch_u,
     ylab="Early Allele - Late Allele",
     xlab="Distance to early allelic max Pause",
     main=paste(organ,"BG", sep=" "),
     las=1, frame=FALSE
)
abline(h=0, col="gray")
for (i in 4:1){
  points(-w:w,target[i,], col = acgt_col[i], type="o", pch=pch_u)
}
legend("topleft", legend=acgt,
       col=acgt_col, bty = "n", lty=1, pch=pch_u)

target=test_early - test_late
organ="3 Tissues" ; target=Tissues3_EarlyPause_1bpapart_KSfdr0.1_early - Tissues3_EarlyPause_1bpapart_KSfdr0.1_late
plot(-w:w, target[1,] , col = acgt_col[1], type="o",
     ylim=c(-0.15,0.15), 
      pch=pch_u,
     ylab="Early Allele - Late Allele",
     xlab="Distance to early allelic max Pause",
     main=paste(organ,"E L at least 1bp apart, KS test fdr<=0.1", sep=" "),
     las=1, frame=FALSE
)
abline(h=0, col="gray")
for (i in 4:1){
  points(-w:w,target[i,], col = acgt_col[i], type="o", pch=pch_u)
}
legend("topleft", legend=acgt,
       col=acgt_col, bty = "n", lty=1, pch=pch_u)

organ="3 Tissues" ; target=Tissues3_LatePause_1bpapart_KSfdr0.1_early - Tissues3_LatePause_1bpapart_KSfdr0.1_late
plot(-w:w, target[1,] , col = acgt_col[1], type="o",
     ylim=c(-0.15,0.15), 
     pch=pch_u,
     ylab="Early Allele - Late Allele",
     xlab="Distance to LATE allelic max Pause",
     main=paste(organ,"E L at least 1bp apart, KS test fdr<=0.1", sep=" "),
     las=1, frame=FALSE
)
abline(h=0, col="gray")
for (i in 4:1){
  points(-w:w,target[i,], col = acgt_col[i], type="o", pch=pch_u)
}
legend("topleft", legend=acgt,
       col=acgt_col, bty = "n", lty=1, pch=pch_u)



