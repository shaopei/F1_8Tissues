setwd("~/Box Sync/KD_IGV/2020July/")

# indel
Tissue="HT"
df<-NULL
for (Tissue in c("HT","SK", "KD")){
df1=read.table(paste(Tissue, "_BothAlleleMaxTSNs_ratio0.5-2_RLmatpat_map3refTomaxTSNmatpat_ClosetIndel.bed", sep = ""))
df1$Tissue=Tissue
df = rbind.data.frame(df,df1)
}
#df= cbind.data.frame(df$Tissue, df[colnames(df)!="Tissue"])

#colnames(df)[7:10]=c("mat_RL", "pat_RL" , "mat_map3RefDist", "pat_map3RefDist" )
View(df)
for (i in 1:dim(df)[1]){
  # pick maxPause
  # if there is a tie, the shorter read length is reported
  df$mat_maxPause_RL[i] = as.numeric(names(sort(table(unlist(strsplit(as.character(df$V7[i]), ","))),decreasing=TRUE)[1]))
  df$pat_maxPause_RL[i] = as.numeric(names(sort(table(unlist(strsplit(as.character(df$V8[i]), ","))),decreasing=TRUE)[1]))
  df$mat_maxPause_map3[i] = as.numeric(names(sort(table(unlist(strsplit(as.character(df$V9[i]), ","))),decreasing=TRUE)[1]))
  df$pat_maxPause_map3[i] = as.numeric(names(sort(table(unlist(strsplit(as.character(df$V10[i]), ","))),decreasing=TRUE)[1]))
  df$mat_Idel_length[i]= length(unlist(strsplit(unlist(strsplit(as.character(df$V14[i]), ","))[1], "")))
  df$pat_Idel_length[i]= length(unlist(strsplit(unlist(strsplit(as.character(df$V14[i]), ","))[2], "")))
  df$maxPauseSite_map3[i] = max(df$mat_maxPause_map3[i],  df$pat_maxPause_map3[i])
  df$maxPauseSite_RL[i] = max(df$mat_maxPause_RL[i],  df$pat_maxPause_RL[i])
}

df$Indel_Len = df$mat_Idel_length - df$pat_Idel_length


for (i in 1:dim(df)[1]){
  m=as.numeric(unlist(strsplit(as.character(df$V7[i]), ",")))
  p=as.numeric(unlist(strsplit(as.character(df$V8[i]), ",")))
  df$RL.p.value[i] = ks.test(m,p) $ p.value
}
for (i in 1:dim(df)[1]){
  m=as.numeric(unlist(strsplit(as.character(df$V9[i]), ",")))
  p=as.numeric(unlist(strsplit(as.character(df$V10[i]), ",")))
  df$map3.p.value[i] = ks.test(m,p) $ p.value
}

df$RL.p.value.fdr = p.adjust(df$RL.p.value, method = "fdr")
df$map3.p.value.fdr = p.adjust(df$map3.p.value, method = "fdr")

# indel length VS map3 position difference
# V15 is the distance of the indel to maxTSN
lim=40


## average
for (i in 1:dim(df)[1]){
  df$mat_AvePause_RL[i] = mean(as.numeric(unlist(strsplit(as.character(df$V7[i]), ","))))
  df$pat_AvePause_RL[i] = mean(as.numeric(unlist(strsplit(as.character(df$V8[i]), ","))))
  df$mat_AvePause_map3[i] = mean(as.numeric(unlist(strsplit(as.character(df$V9[i]), ","))))
  df$pat_AvePause_map3[i] = mean(as.numeric(unlist(strsplit(as.character(df$V10[i]), ","))))
}

df$target=df$V15 <= df$maxPauseSite_map3 & (df$map3.p.value.fdr <=0.1 )
dim(df)
sum(df$target)
sum(df$map3.p.value.fdr<0.1)
sum(df$map3.p.value.fdr>0.1)
# plot( (df$mat_AvePause_map3- df$pat_AvePause_map3)[df$target], 
#       (df$mat_Idel_length - df$pat_Idel_length)[df$target],
#       xlim=c(-1*lim,lim), ylim=c(-1*lim,lim),
#       xlab="B6 - CAST pause position",
#       ylab="B6 - CAST indel length",
#       frame=F, 
#       pch=19, col=rgb(0,0,0,alpha = 0.25))
pch_u=19
plot( (df$mat_AvePause_map3- df$pat_AvePause_map3)[df$target& df$Tissue=="HT"],
        (df$mat_Idel_length - df$pat_Idel_length)[df$target & df$Tissue=="HT"],
        xlim=c(-1*lim,lim),
      ylim=c(-1*lim,lim),
      xlab="B6 - CAST pause position",
      ylab="B6 - CAST indel length",
      pch=pch_u, frame=F,
      col="red")
# points( (df$mat_AvePause_map3- df$pat_AvePause_map3)[df$target& df$Tissue=="HT"], 
#         (df$mat_Idel_length - df$pat_Idel_length)[df$target & df$Tissue=="HT"],
#         pch=pch_u, 
#         col="red")
points( (df$mat_AvePause_map3- df$pat_AvePause_map3)[df$target& df$Tissue=="SK"], 
        (df$mat_Idel_length - df$pat_Idel_length)[df$target & df$Tissue=="SK"],
        pch=pch_u, 
        col="dark orange")
points( (df$mat_AvePause_map3- df$pat_AvePause_map3)[df$target& df$Tissue=="KD"], 
      (df$mat_Idel_length - df$pat_Idel_length)[df$target & df$Tissue=="KD"],
      pch=pch_u, col="blue")

abline(a=5, b=-1, col="gray")
abline(a=0,b=-1, col="dark gray")
abline(a=-5, b=-1, col="gray")

legend("topright", 
       legend = c("HT", "SK", "KD"),
       pch=pch_u,
#       cex=2, 
       lty=c(0,0),
       #bty="n",
       #lwd=1.5, 
       #density=c(10000,25),
       #angle=c(180,45),
       #angle=45,
       col=c("red", "dark orange", "blue")
       , bty = "n"
)


# get read length distribution
for (Tissue in c("HT","SK", "KD")){
  pdf(file = paste(Tissue, "_read_length_distribution.pdf"))
ide_name=paste(Tissue,"_mat_identical_read_length.txt",sep="")
mat_name=paste(Tissue,"_mat_specific_read_length.txt",sep="")
pat_name=paste(Tissue,"_pat_specific_read_length.txt",sep="")

ide<-read.table(ide_name)
mat<-read.table(mat_name)
pat<-read.table(pat_name)

par(mfrow=c(3,1))
xlim_top=101
hist(mat$V1[mat$V1<=100], 
     xlim = c(0, xlim_top), las=1, 
     breaks = seq(0.5, xlim_top,1),
     main=Tissue)
hist(pat$V1[pat$V1<=100], 
     breaks = seq(0.5, xlim_top,1),
     xlim = c(0, xlim_top), las=1)
hist(ide$V1[ide$V1<=100], 
     breaks = seq(0.5, xlim_top,1),
     xlim = c(0, xlim_top), las=1)
dev.off()
}
