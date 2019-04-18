#R --vanilla --slave --args $(pwd) Tissue < Find_consistent_blocks.R

#arguments here
args=(commandArgs(TRUE))
setwd(args[1])
Tissue=args[2]

#setwd("~/Box Sync/Danko_lab_work/F1_8Tissues/consistent_block")
library("metap")

write.bed<-function ( df.bed, file.bed, compress=FALSE ) {
  #options("scipen"=100, "digits"=4);
  temp <- tempfile(fileext=".bed");
  write.table( df.bed, file=temp, quote=F, row.names=F, col.names=F, sep="\t");
  if(compress)
    system(paste0("sort-bed ", temp,  " | bgzip > ",  file.bed ))
  else
    system(paste0("sort-bed ", temp,  " > ",  file.bed ));
  invisible(unlink(temp));
}


MB6_sumlog <- function(Tissue, strand){
  Group=paste(Tissue,"MB6", sep="_")
  df=read.table(paste(Group,"all_R1_gencode.vM20.annotation_transcript.merged_cov_binomtest.bed", sep = "_"))
  colnames(df)[grep("V4",colnames(df))]="win_count_all_samples"
  colnames(df)[grep("V9",colnames(df))]="p_value_all_samples"
  #colnames(df)[grep("V12",colnames(df))]="V12.all"
  col_to_keep=c("V1","V2", "V3", "V10", "V11", "V12","win_count_all_samples","p_value_all_samples")
  col_to_merge=c("V1","V2", "V3", "V10", "V11", "V12")
  df=df[,col_to_keep]
  #df=df[,1:3]
  for (pool in c("A", "F", "G")){
    sub_f=read.table(file = paste(Group, pool, "R1_gencode.vM20.annotation_transcript.merged_cov_binomtest.bed", sep = "_"))
    colnames(sub_f)[grep("V4",colnames(sub_f))]=paste("V4", pool, sep=".")
    colnames(sub_f)[grep("V9",colnames(sub_f))]=paste("V9", pool, sep=".")
    #colnames(sub_f)[grep("V12",colnames(sub_f))]=paste("V12", pool, sep=".")
    df=merge(df,sub_f, by =col_to_merge, suffixes =c(".x",paste(".",pool, sep="")))
    # only keep the rows that are present in all samples
    df=df[,c(col_to_keep,colnames(df)[grep("V4",colnames(df))],colnames(df)[grep("V9",colnames(df))])]
  }
  df[,grep("V9",colnames(df))][df[,grep("V9",colnames(df))]==0]=1e-20
  result = unlist( apply( df[,grep("V9",colnames(df))],1, function(x) { r<-sumlog(x); r$p;}) )
  df$sumlog=result
  write.bed(df[df$sumlog<=0.05,], file=paste(Group, pool, "R1_gencode.vM20.transcript_ABconsistent_FisherMethodP0.05.bed", sep = "_"))
  write.bed(df, file=paste(Group, pool, "R1_gencode.vM20.transcript_ABconsistent.bed", sep = "_"))
}



PB6_sumlog <- function(Tissue, strand){
  Group=paste(Tissue,"PB6", sep="_")
  df=read.table(paste(Group,"all_R1_HMM", strand, "agreeCount.bed", sep = "_"))
  df=df[,1:3]
  for (pool in c("B", "C", "D", "E")){
    sub_f=read.table(file = paste(Group, pool, "R1_HMM",strand, "binomtest.bed", sep = "_"))
    df=merge(df,sub_f, by =c("V1","V2", "V3"),suffixes =c(".x",paste(".",pool, sep="")))
  }
  df=df[,c("V1","V2", "V3",colnames(df)[grep("V9",colnames(df))])]
  df[,4:dim(df)[2]][df[,4:dim(df)[2]]==0]=1e-20
  result = unlist( apply( df[,4:dim(df)[2]],1, function(x) { r<-sumlog(x); r$p;}) )
  df$sumlog=result
  write.bed(df[df$sumlog<=0.05,], file=paste(Group,"all_R1_HMM", strand, "agreeCount_strict_P0.05.bed", sep = "_"))
  write.bed(df, file=paste(Group,"all_R1_HMM", strand, "agreeCount_strict.bed", sep = "_"))
}


MB6_sumlog(Tissue, "plus")
MB6_sumlog(Tissue, "minus")
PB6_sumlog(Tissue, "plus")
PB6_sumlog(Tissue, "minus")
  
