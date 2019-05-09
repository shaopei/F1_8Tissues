library(UpSetR)
setwd("~/Box Sync/Danko_lab_work/F1_8Tissues/UpSetR")
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



# imprinting cluster
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


