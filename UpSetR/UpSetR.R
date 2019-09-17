library(UpSetR)
setwd("~/Box Sync/Danko_lab_work/F1_8Tissues/UpSetR")
par(mar=c(6.1, 7.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
par(mgp=c(3,1,0))
par(cex.lab=2.2, cex.axis=2.2)
# imprinting cluster regardless of strandness
Tissue_list=c( "BN","SP","HT","SK","KD","ST","GI","LV")
df=read.table("T8_2Strand_p0.05_effect_imprinting.bed_cluster", header = F)
#df=read.table("T8_2Strand_1FMp0.05_effect_imprinting.bed_cluster", header = F)
#df=read.table("T8_2Strand_p0.1_effect_imprinting.bed_cluster", header = F)
for (kkk in Tissue_list){
  df$tmp=0
  df$tmp[grepl(kkk, df$V6)]=1
  colnames(df)[grep("tmp", colnames(df))]=kkk
}
upset(df, nsets = 8, sets =Tissue_list, keep.order = T, order.by = "degree")
df$TissueCounts= rowSums(df[ , match(Tissue_list , names(df) ) ] )  
a=hist(df$TissueCounts[df$TissueCounts>1], breaks = seq(-0.5,9,1))
a=hist(df$TissueCounts[df$TissueCounts >=1], breaks = seq(0.5,9,1),
       xlab="Number of organs", ylab="Number of clusters", main="",
       las=1, freq = T)
a$counts

hist(df$TissueCounts, breaks = seq(-0.5,9,1), col="blue",las=1, ylim=c(0,80))
hist(df$TissueCounts, breaks = seq(-0.5,9,1), col="blue",las=1,add=T)
hist(df$TissueCounts, breaks = seq(-0.5,9,1), col="red", density=25,angle=45,add=T)
hist(df$TissueCounts, breaks = seq(-0.5,9,1), col="red", density=25,angle=45)
sum(df$TissueCounts==1)

legend("topright", 
       legend = c("100%", "25%"),
       #pch=c(15,15),
       cex=2, 
       lty=c(0,0),
       #bty="n",
       lwd=1.5, 
       density=c(25,10000),
       angle=c(45,180),
       #angle=45,
       fill=c("blue", "yellow")
       , bty = "n"
)

# strain effect cluster regardless of strandness
Tissue_list=c( "BN","SP","HT","SK","KD","ST","GI","LV")
df=read.table("T8_2Strand_p0.05_effect_strain.bed_cluster", header = F)
for (kkk in Tissue_list){
  df$tmp=0
  df$tmp[grepl(kkk, df$V6)]=1
  colnames(df)[grep("tmp", colnames(df))]=kkk
}
upset(df, nsets = 8, sets =Tissue_list, keep.order = T, order.by = "degree", nintersects=100)
df$TissueCounts= rowSums(df[ , match(Tissue_list , names(df) ) ] )  
a=hist(df$TissueCounts, breaks = seq(-0.5,9,1))
a=hist(df$TissueCounts[df$TissueCounts >=1], breaks = seq(0.5,9,1),
       xlab="Number of organs", ylab="Number of clusters", main="",
       las=1, freq = T)
barplot(a$counts, names.arg=seq(1,8,1))

hist(df$TissueCounts, breaks = seq(-0.5,9,1), col="blue",las=1, ylim=c(0,80))
hist(df$TissueCounts, breaks = seq(-0.5,9,1), col="blue",las=1,add=T)
hist(df$TissueCounts, breaks = seq(-0.5,9,1), col="red", density=25,angle=45,add=T)
hist(df$TissueCounts, breaks = seq(-0.5,9,1), col="red", density=25,angle=45)
sum(df$TissueCounts==1)

legend("topright", 
       legend = c("100%", "25%"),
       #pch=c(15,15),
       cex=2, 
       lty=c(0,0),
       #bty="n",
       lwd=1.5, 
       density=c(25,10000),
       angle=c(45,180),
       #angle=45,
       fill=c("blue", "yellow")
       , bty = "n"
)

# cluster combined strain effect and imprinting regardless of strandness
Tissue_list=c( "BN","SP","HT","SK","KD","ST","GI","LV")
df=read.table("T8_2Strand_p0.05_effect_strain+imprinting.bed_cluster", header = F)
for (kkk in Tissue_list){
  df$tmp=0
  df$tmp[grepl(kkk, df$V6)]=1
  colnames(df)[grep("tmp", colnames(df))]=kkk
}
upset(df, nsets = 8, sets =Tissue_list, keep.order = T, order.by = "degree", nintersects=100)
df$TissueCounts= rowSums(df[ , match(Tissue_list , names(df) ) ] )  
a=hist(df$TissueCounts[df$TissueCounts >=1], breaks = seq(0.5,9,1),
       xlab="Number of organs", ylab="Number of clusters", main="",
       las=1, freq = T)
barplot(a$counts, names.arg=seq(1,8,1))

hist(df$TissueCounts, breaks = seq(-0.5,9,1), col="blue",las=1, ylim=c(0,80))
hist(df$TissueCounts, breaks = seq(-0.5,9,1), col="blue",las=1,add=T)
hist(df$TissueCounts, breaks = seq(-0.5,9,1), col="red", density=25,angle=45,add=T)
hist(df$TissueCounts, breaks = seq(-0.5,9,1), col="red", density=25,angle=45)
sum(df$TissueCounts==1)



# imprinting genes/ncRNA merged regardless of strandness
Tissue_list=c( "BN","SP","HT","SK","KD","ST","GI","LV")
df=read.table("T8_ncRNA_p0.05_effect_strain_merged.bed", header = F)
for (kkk in Tissue_list){
  df$tmp=0
  df$tmp[grepl(kkk, df$V4)]=1
  colnames(df)[grep("tmp", colnames(df))]=kkk
}
upset(df, nsets = 8, sets =Tissue_list, keep.order = T, order.by = "degree")
df$TissueCounts= rowSums(df[ , match(Tissue_list , names(df) ) ] )  

a=hist(df$TissueCounts[df$TissueCounts >=1], breaks = seq(0.5,9,1),
       xlab="Number of organs", ylab="Number of clusters", main="",
       las=1, freq = T)
a$counts

hist(df$TissueCounts, breaks = seq(-0.5,9,1), col="pink",las=1, main="ncRNA, strain")

df=read.table("T8_ProteinCodingGene_p0.05_effect_strain_merged.bed", header = F)
for (kkk in Tissue_list){
  df$tmp=0
  df$tmp[grepl(kkk, df$V4)]=1
  colnames(df)[grep("tmp", colnames(df))]=kkk
}
upset(df, nsets = 8, sets =Tissue_list, keep.order = T, order.by = "degree")
df$TissueCounts= rowSums(df[ , match(Tissue_list , names(df) ) ] )  

a=hist(df$TissueCounts[df$TissueCounts >=1], breaks = seq(0.5,9,1),
       xlab="Number of organs", ylab="Number of clusters", main="",
       las=1, freq = T)
a$counts

hist(df$TissueCounts, breaks = seq(-0.5,9,1), col="red",las=1)




hist(df$TissueCounts, breaks = seq(-0.5,9,1), col="blue",las=1,add=T)
hist(df$TissueCounts, breaks = seq(-0.5,9,1), col="red", density=25,angle=45,add=T)
hist(df$TissueCounts, breaks = seq(-0.5,9,1), col="red", density=25,angle=45)
sum(df$TissueCounts==1)

legend("topright", 
       legend = c("100%", "25%"),
       #pch=c(15,15),
       cex=2, 
       lty=c(0,0),
       #bty="n",
       lwd=1.5, 
       density=c(25,10000),
       angle=c(45,180),
       #angle=45,
       fill=c("blue", "yellow")
       , bty = "n"
)




























##### not use #####
# imprnting blocks
mdf=read.table("T8_Merge_minus_p0.05_effect_imprinting.bed", header = T,  comment.char = "/")
mdf$name=paste(mdf$X.chrm,mdf$chrStart,mdf$chrEnd,"-",sep="_")
pdf=read.table("T8_Merge_plus_p0.05_effect_imprinting.bed", header = T,  comment.char = "/")
pdf$name=paste(pdf$X.chrm,pdf$chrStart,pdf$chrEnd,"+",sep="_")
df=rbind(pdf,mdf)
df=df[,c(grep("name",colnames(df)),grep("MB6",colnames(df)))]
df[,grep("B6",colnames(df))]=ifelse(df[,grep("B6",colnames(df))]=="EMPTY", 0, 1)
colnames(df)= unlist(lapply(colnames(df), function(x) { r<-strsplit(x, "_"); unlist(r)[1]}))
upset(df, nsets = 8, sets =c( "BN","HT","SK","SP","LG","LV","GI","ST"), keep.order = T, order.by = "degree")


# imprinting cluster
mdf=read.table("T8_minus_p0.05_effect_imprinting.bed_cluster", header = F)
mdf$name=paste(mdf$V1,mdf$V2,mdf$V3,"-",sep="_")
pdf=read.table("T8_plus_p0.05_effect_imprinting.bed_cluster", header = F)
pdf$name=paste(pdf$V1,pdf$V2,pdf$V3,"+",sep="_")
df=rbind(pdf,mdf)
for (kkk in c( "BN","HT","SK","SP","LG","LV","GI","ST")){
  df$tmp=0
  df$tmp[grepl(kkk, df$V4)]=1
  colnames(df)[grep("tmp", colnames(df))]=kkk
}
upset(df, nsets = 8, sets =c( "BN","HT","SK","SP","LG","LV","GI","ST"), keep.order = T, order.by = "degree")



# strain cluster
mdf=read.table("T8_minus_p0.05_effect_strain.bed_cluster", header = F)
mdf$name=paste(mdf$V1,mdf$V2,mdf$V3,"-",sep="_")
pdf=read.table("T8_plus_p0.05_effect_strain.bed_cluster", header = F)
pdf$name=paste(pdf$V1,pdf$V2,pdf$V3,"+",sep="_")
df=rbind(pdf,mdf)
for (kkk in c( "BN","HT","SK","SP","LG","LV","GI","ST")){
  df$tmp=0
  df$tmp[grepl(kkk, df$V4)]=1
  colnames(df)[grep("tmp", colnames(df))]=kkk
}
upset(df, nsets = 8, sets =c( "BN","HT","SK","SP","LG","LV","GI","ST"), keep.order = T, order.by = "degree")



# imprinting genes /ncRNA