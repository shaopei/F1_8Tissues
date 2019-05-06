library(UpSetR)
mdf=read.table("T8_Merge_minus_p0.05_effect_imprinting.bed", header = T,  comment.char = "/")
mdf$name=paste(mdf$X.chrm,mdf$chrStart,mdf$chrEnd,"-",sep="_")
pdf=read.table("T8_Merge_plus_p0.05_effect_imprinting.bed", header = T,  comment.char = "/")
pdf$name=paste(pdf$X.chrm,pdf$chrStart,pdf$chrEnd,"+",sep="_")
df=rbind(pdf,mdf)
df=df[,c(grep("name",colnames(df)),grep("MB6",colnames(df)))]
df[,grep("B6",colnames(df))]=ifelse(df[,grep("B6",colnames(df))]=="EMPTY", 0, 1)
colnames(df)= unlist(lapply(colnames(df), function(x) { r<-strsplit(x, "_"); unlist(r)[1]}))
upset(df, nsets = 8, sets =c( "BN","HT","SK","SP","LG","LV","GI","ST"), keep.order = T, order.by = "degree")


