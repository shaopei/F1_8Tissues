#R --vanilla --slave --args $(pwd) Tissue < KStest.R

#arguments here
args=(commandArgs(TRUE))
setwd(args[1])
Tissue=args[2]
# if (length(args) <3){
#   pvlaue_cutoff=0.05
# }else{
#   pvlaue_cutoff=as.numeric(args[3])
# }

#setwd("~/Box Sync/Danko_lab_work/F1_8Tissues/Kidney_and_SingleRunOn")
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


mat<-read.table(paste(Tissue, "gencode.vM20.annotation_transcript_30bp_wSNP_10+reads_mat.perBase.bed", sep = "_"))
pat<-read.table(paste(Tissue, "gencode.vM20.annotation_transcript_30bp_wSNP_10+reads_pat.perBase.bed", sep = "_"))
name <- read.table(paste(Tissue, "gencode.vM20.annotation_transcript_30bp_wSNP_10+reads_uniq.bed", sep = "_"))

for (i in 1:dim(mat)[1]){
name$p.value[i] = ks.test(as.numeric(mat[i,]),as.numeric(pat[i,])) $ p.value
}
name$p.value.fdr = p.adjust(name$p.value, method = "fdr")

name=name[order(name$p.value, decreasing = FALSE),]

write.table(name, file=paste(Tissue, "gencode.vM20.annotation_transcript_30bp_wSNP_10+reads_uniq_pValue.bed", sep = "_"), quote=F, row.names=F, col.names=F, sep="\t")
#write.bed(name, paste(Tissue, "gencode.vM20.annotation_transcript_30bp_wSNP_10+reads_uniq_pValue.bed", sep = "_"))

