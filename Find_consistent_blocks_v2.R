#R --vanilla --slave --args $(pwd) Tissue < Find_consistent_blocks_v2.R

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
  # use BN_MB6_all_R1_HMM_plus_binomtest.bed for df
  df=read.table(paste(Group,"all_R1_HMM", strand, "binomtest.bed", sep = "_"))
  df$winP = substr(df$V4,1,1)
  colnames(df)[grep("V4",colnames(df))]="win_count_all_samples"
  colnames(df)[grep("V9",colnames(df))]="p_value_all_samples"
  col_to_keep=c("V1","V2", "V3", "win_count_all_samples","p_value_all_samples","winP")
  col_to_merge=c("V1","V2", "V3", "winP")
  
  df=df[,col_to_keep]
  #df=df[,1:3]
  for (pool in c("A", "F", "G")){
    #BN_MB6_G_R1_HMM_plus_binomtest.bed
    file = paste(Group, pool, "R1_HMM",strand, "binomtest.bed", sep = "_")
    if (file.exists(file)){
      sub_f=read.table(file = file)
      sub_f$winP = substr(sub_f$V4,1,1)
      #colnames(sub_f)[grep("V4",colnames(sub_f))]=paste("V4", pool, sep=".")
      colnames(sub_f)[grep("V9",colnames(sub_f))]=paste("V9", pool, sep=".")
      #colnames(sub_f)[grep("V12",colnames(sub_f))]=paste("V12", pool, sep=".")
      df=merge(df,sub_f, by =col_to_merge, suffixes =c(".x",paste(".",pool, sep="")))
      # only keep the rows that are present in all samples
      df=df[,c(col_to_keep,colnames(df)[grep("V9",colnames(df))])]
    }
  }
  df[,grep("V9",colnames(df))][df[,grep("V9",colnames(df))]==0]=1e-20
  result = unlist( apply( df[,grep("V9",colnames(df))],1, function(x) { r<-sumlog(x); r$p;}) )
  df$sumlog=result
  df=df[,c(col_to_keep,"sumlog", colnames(df)[grep("V9",colnames(df))])]
  write.bed(df[df$sumlog<=0.05,], file=paste(Group,"combined_R1_HMM", strand, "ABconsistent_FisherMethodP0.05.bed", sep = "_"))
  write.bed(df, file=paste(Group,"combined_R1_HMM", strand, "ABconsistent.bed", sep = "_"))
}
  


PB6_sumlog <- function(Tissue, strand){
  Group=paste(Tissue,"PB6", sep="_")
  # use BN_PB6_all_R1_HMM_plus_binomtest.bed for df
  df=read.table(paste(Group,"all_R1_HMM", strand, "binomtest.bed", sep = "_"))
  df$winP = substr(df$V4,1,1)
  colnames(df)[grep("V4",colnames(df))]="win_count_all_samples"
  colnames(df)[grep("V9",colnames(df))]="p_value_all_samples"
  col_to_keep=c("V1","V2", "V3", "win_count_all_samples","p_value_all_samples","winP")
  col_to_merge=c("V1","V2", "V3", "winP")
  
  df=df[,col_to_keep]
  #df=df[,1:3]
  for (pool in c("B", "C", "D", "E")){
    #BN_PB6_B_R1_HMM_plus_binomtest.bed
    file = paste(Group, pool, "R1_HMM",strand, "binomtest.bed", sep = "_")
    if (file.exists(file)){
      sub_f=read.table(file = file)
      sub_f$winP = substr(sub_f$V4,1,1)
      #colnames(sub_f)[grep("V4",colnames(sub_f))]=paste("V4", pool, sep=".")
      colnames(sub_f)[grep("V9",colnames(sub_f))]=paste("V9", pool, sep=".")
      #colnames(sub_f)[grep("V12",colnames(sub_f))]=paste("V12", pool, sep=".")
      df=merge(df,sub_f, by =col_to_merge, suffixes =c(".x",paste(".",pool, sep="")))
      # only keep the rows that are present in all samples
      df=df[,c(col_to_keep,colnames(df)[grep("V9",colnames(df))])]
    }
  }
  df[,grep("V9",colnames(df))][df[,grep("V9",colnames(df))]==0]=1e-20
  result = unlist( apply( df[,grep("V9",colnames(df))],1, function(x) { r<-sumlog(x); r$p;}) )
  df$sumlog=result
  df=df[,c(col_to_keep,"sumlog", colnames(df)[grep("V9",colnames(df))])]
  write.bed(df[df$sumlog<=0.05,], file=paste(Group,"combined_R1_HMM", strand, "ABconsistent_FisherMethodP0.05.bed", sep = "_"))
  write.bed(df, file=paste(Group,"combined_R1_HMM", strand, "ABconsistent.bed", sep = "_"))
}


MB6_sumlog(Tissue, "plus")
MB6_sumlog(Tissue, "minus")
PB6_sumlog(Tissue, "plus")
PB6_sumlog(Tissue, "minus")
  
