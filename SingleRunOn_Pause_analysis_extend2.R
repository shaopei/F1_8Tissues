setwd("~/Box Sync/KD_IGV/2020July/")
df = read.table("HT_matReads_patReads_TSS_maxTSNs_ratio0.5-2_ReadLengthPValue_map3TomaxTSNPValue.bed")
colnames(df)[8] = "ReadLenKSfdr"
colnames(df)[11] = "Map3map2refDistMaxTSNKSfdr"
plot(df$ReadLenKSfdr, df$Map3map2refDistMaxTSNKSfdr,
     pch=19, col=rgb(0,0,0,alpha = 0.125))


dim(df)
sum(df$ReadLenKSfdr<0.1 & df$Map3map2refDistMaxTSNKSfdr>0.9)
sum(df$ReadLenKSfdr>0.9 & df$Map3map2refDistMaxTSNKSfdr<0.1)
sum(df$ReadLenKSfdr<0.1 & df$Map3map2refDistMaxTSNKSfdr<0.1)


ReadLenKS_df=df[df$ReadLenKSfdr<0.1,]
View(ReadLenKS_df)

Map3Distance_df=df[df$Map3map2refDistMaxTSNKSfdr<0.1,]
View(Map3Distance_df)


mat<-read.table("HT_matReads_patReads_TSS_maxTSNs_ratio0.5-2_mat.ReadLength.bed")
pat<-read.table("HT_matReads_patReads_TSS_maxTSNs_ratio0.5-2_pat.ReadLength.bed")

mat_readLength <- NULL
pat_readLength <- NULL
for (i in 1:dim(mat)[1]){
  mat_readLength = c( mat_readLength, as.numeric(strsplit(as.character(mat$V7[i]), ",")[[1]]))
  pat_readLength = c( pat_readLength, as.numeric(strsplit(as.character(pat$V7[i]), ",")[[1]]))
}

save.image("x.RData")
setwd("~/Box Sync/KD_IGV/2020July/")
load("x.RData")
readLength=c(mat_readLength, pat_readLength)
hist(readLength, breaks = seq(-0.5,100,1))
hist(mat_readLength, breaks = seq(-0.5,100,1))
hist(pat_readLength, breaks = seq(-0.5,100,1))


# if there is a tie, the shorter read length is reported
for (i in 1:dim(mat)[1]){
  mat$maxPause_RL[i] = as.numeric(names(sort(table(unlist(strsplit(as.character(mat$V7[i]), ","))),decreasing=TRUE)[1]))
  pat$maxPause_RL[i] = as.numeric(names(sort(table(unlist(strsplit(as.character(pat$V7[i]), ","))),decreasing=TRUE)[1]))
}

combine = mat
colnames(combine)[colnames(combine)=="maxPause_RL"] = "mat_maxPause_RL"
combine$pat_maxPause_RL = pat$maxPause_RL[pat$V2 == combine$V2]

for (i in 1:dim(mat)[1]){
  m=as.numeric(strsplit(as.character(mat$V7[i]), ",")[[1]])
  p=as.numeric(strsplit(as.character(pat$V7[i]), ",")[[1]])
  combine$p.value[i] = ks.test(m,p) $ p.value
}
combine$p.value.fdr = p.adjust(combine$p.value, method = "fdr")
combine$deltaDistMaxPause_RL = combine$mat_maxPause_RL - combine$pat_maxPause_RL

View(combine)

hist(abs(combine$deltaDistMaxPause)[combine$p.value.fdr > 0.9] # & combine$deltaDistMaxPause != 0], 
     ,breaks = seq(-0.5,30,1),
     freq = F,
     col = "blue",
     xlab="Difference between allelic read length (highest frequency) from same maxTSN",
     main="Read Length difference"
)
hist(abs(combine$deltaDistMaxPause)[combine$p.value.fdr <= 0.1] #  & combine$deltaDistMaxPause != 0], 
     ,breaks = seq(-0.5,30,1),
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


df=read.table("HT_BothAlleleMaxTSNs_ratio0.5-2_RLmatpat_map3refTomaxTSNmatpat.bed")
#colnames(df)[7:10]=c("mat_RL", "pat_RL" , "mat_map3RefDist", "pat_map3RefDist" )
View(df)
for (i in 1:dim(df)[1]){
  # pick maxPause
  # if there is a tie, the shorter read length is reported
  df$mat_maxPause_RL[i] = as.numeric(names(sort(table(unlist(strsplit(as.character(df$V7[i]), ","))),decreasing=TRUE)[1]))
  df$pat_maxPause_RL[i] = as.numeric(names(sort(table(unlist(strsplit(as.character(df$V8[i]), ","))),decreasing=TRUE)[1]))
  df$mat_maxPause_map3[i] = as.numeric(names(sort(table(unlist(strsplit(as.character(df$V9[i]), ","))),decreasing=TRUE)[1]))
  df$pat_maxPause_map3[i] = as.numeric(names(sort(table(unlist(strsplit(as.character(df$V10[i]), ","))),decreasing=TRUE)[1]))
    }
plot(df$mat_maxPause_RL- df$pat_maxPause_RL, df$mat_maxPause_map3- df$pat_maxPause_map3,
     xlim=c(-30,30), ylim=c(-30,30),
     pch=19, col=rgb(0,0,0,alpha = 0.125))


plot(df$mat_maxPause_RL- df$pat_maxPause_RL, df$mat_maxPause_map3- df$pat_maxPause_map3,
     xlim=c(-30,30), ylim=c(-30,30),
     pch=19, col=rgb(0,0,0,alpha = 0.125))
df$pick = df$RL.p.value.fdr <0.1
points((df$mat_maxPause_RL- df$pat_maxPause_RL)[df$pick]  , (df$mat_maxPause_map3- df$pat_maxPause_map3)[df$pick],
       col="red"
       #xlim=c(-20,20), ylim=c(-20,20)
)
df$pick = df$map3.p.value.fdr <0.1
points((df$mat_maxPause_RL- df$pat_maxPause_RL)[df$pick]  , (df$mat_maxPause_map3- df$pat_maxPause_map3)[df$pick],
       col="blue"
       #xlim=c(-20,20), ylim=c(-20,20)
)
df$pick = df$RL.p.value.fdr <0.1 & df$map3.p.value.fdr <0.1
points((df$mat_maxPause_RL- df$pat_maxPause_RL)[df$pick]  , (df$mat_maxPause_map3- df$pat_maxPause_map3)[df$pick],
       col="purple"
       #xlim=c(-20,20), ylim=c(-20,20)
)
abline(v=0)
abline(h=0)

for (i in 1:dim(df)[1]){
  df$mat_AvePause_RL[i] = mean(as.numeric(unlist(strsplit(as.character(df$V7[i]), ","))))
  df$pat_AvePause_RL[i] = mean(as.numeric(unlist(strsplit(as.character(df$V8[i]), ","))))
  df$mat_AvePause_map3[i] = mean(as.numeric(unlist(strsplit(as.character(df$V9[i]), ","))))
  df$pat_AvePause_map3[i] = mean(as.numeric(unlist(strsplit(as.character(df$V10[i]), ","))))
}
plot((df$mat_AvePause_RL - df$pat_AvePause_RL) , (df$mat_AvePause_map3- df$pat_AvePause_map3),
     xlim=c(-20,20), ylim=c(-20,20),
     pch=19, col=rgb(0,0,0,alpha = 0.125))
#abline(v=0)
#abline(h=0)

df$pick = df$RL.p.value.fdr <0.1
points(abs(df$mat_AvePause_RL - df$pat_AvePause_RL)[df$pick]  , abs(df$mat_AvePause_map3- df$pat_AvePause_map3)[df$pick],
       col="red"
       #xlim=c(-20,20), ylim=c(-20,20)
)
df$pick = df$map3.p.value.fdr <0.1
points(abs(df$mat_AvePause_RL - df$pat_AvePause_RL)[df$pick]  , abs(df$mat_AvePause_map3- df$pat_AvePause_map3)[df$pick],
       col="blue"
       #xlim=c(-20,20), ylim=c(-20,20)
)

df$pick = df$map3.p.value.fdr <0.1 & df$RL.p.value.fdr <0.1
points((df$mat_AvePause_RL - df$pat_AvePause_RL)[df$pick]  , (df$mat_AvePause_map3- df$pat_AvePause_map3)[df$pick],
       col="purple"
       #xlim=c(-20,20), ylim=c(-20,20)
)
abline(v=0)
abline(h=0)
abline(a=0, b=1)
abline(a=1, b=1)



 # KS test result?
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
plot(df$RL.p.value.fdr, df$map3.p.value.fdr,
     pch=19, col=rgb(0,0,0,alpha = 0.125))



u=max(abs(df$mat_maxPause_RL - df$pat_maxPause_RL))+1
hist(abs(df$mat_maxPause_RL - df$pat_maxPause_RL) # & combine$deltaDistMaxPause != 0], 
     ,breaks = seq(-0.5,u,1),
     freq = F,
     col = "blue",
     xlab="Difference between allelic read length (highest frequency) from same maxTSN",
     main= paste(Tissue, "Read Length difference", sep=" ")
)
hist(abs(df$mat_maxPause_RL - df$pat_maxPause_RL)[df$RL.p.value.fdr<=0.1] #  & combine$deltaDistMaxPause != 0], 
     ,breaks = seq(-0.5,u,1),
     freq = F,
     density = 50, col="red",
     add=T
)

u=max(abs(df$mat_maxPause_map3 - df$pat_maxPause_map3))+1
hist(abs(df$mat_maxPause_map3 - df$pat_maxPause_map3) # & combine$deltaDistMaxPause != 0], 
     ,breaks = seq(-0.5,u,1),
     freq = F,
     col = "blue",
     xlab="Difference between allelic read length (highest frequency) from same maxTSN",
     main= paste(Tissue, "Read Length difference", sep=" ")
)
hist(abs(df$mat_maxPause_map3 - df$pat_maxPause_map3)[df$map3.p.value.fdr <=0.1] #  & combine$deltaDistMaxPause != 0], 
     ,breaks = seq(-0.5,u,1),
     freq = F,
     density = 50, col="red",
     add=T
)


u=max(abs(df$mat_maxPause_map3 - df$pat_maxPause_map3 - (df$mat_maxPause_RL - df$pat_maxPause_RL) ))+1
hist(abs(df$mat_maxPause_map3 - df$pat_maxPause_map3 - (df$mat_maxPause_RL - df$pat_maxPause_RL) )
     ,breaks = seq(-0.5,u,1),
     freq = F,
     col = "blue",
)
hist(abs(df$mat_maxPause_map3 - df$pat_maxPause_map3 - (df$mat_maxPause_RL - df$pat_maxPause_RL))[df$map3.p.value.fdr <=0.1 | df$RL.p.value.fdr <=0.1] #  & combine$deltaDistMaxPause != 0], 
     ,breaks = seq(-0.5,u,1),
     freq = F,
     density = 50, col="red",
     add=T
)
legend("topright", 
       legend = c("All", "RL or map3TomaxTSNs KS FDR <= 0.1"),
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

u=max(abs(df$mat_AvePause_map3 - df$pat_AvePause_map3 - (df$mat_AvePause_RL - df$pat_AvePause_RL) ))+1
hist(abs(df$mat_AvePause_map3 - df$pat_AvePause_map3 - (df$mat_AvePause_RL - df$pat_AvePause_RL) )
     ,breaks = seq(-0.5,u,1), 
     ylim=c(0,10),
     #freq = F,
     col = "blue",
)
hist(abs(df$mat_AvePause_map3 - df$pat_AvePause_map3 - (df$mat_AvePause_RL - df$pat_AvePause_RL))[df$map3.p.value.fdr <=0.1 | df$RL.p.value.fdr <=0.1] 
     ,breaks = seq(-0.5,u,1),
     #freq = F,
     density = 50, col="red",
     add=T
)
legend("topright", 
       legend = c("All", "RL or map3TomaxTSNs KS FDR <= 0.1"),
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
