#R --vanilla --slave --args $(pwd) Tissue gencode.vM20.annotation_transcript_100bp_5mat5pat_uniq < KStest.R
# only use bed regions with at least 5 VALID mat AND 5 pat VALID reads

#arguments here
args=(commandArgs(TRUE))
setwd(args[1])
Tissue=args[2]
name_body=args[3]
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


mat<-read.table(paste(Tissue, name_body, "mat.perBase.bed", sep = "_"))
pat<-read.table(paste(Tissue, name_body, "pat.perBase.bed", sep = "_"))
name <- read.table(paste(paste(Tissue, name_body, sep = "_"), "bed", sep ="."))

for (i in 1:dim(mat)[1]){
  mat$read.count[i] = length(as.numeric(strsplit(as.character(mat$V7[i]), ",")[[1]]))
  pat$read.count[i] = length(as.numeric(strsplit(as.character(pat$V7[i]), ",")[[1]]))
  }

new_mat = mat[(mat$read.count>=5 & pat$read.count>=5),]
new_pat = pat[(mat$read.count>=5 & pat$read.count>=5),]
name = name[(mat$read.count>=5 & pat$read.count>=5),]

for (i in 1:dim(new_mat)[1]){
  m=as.numeric(strsplit(as.character(new_mat$V7[i]), ",")[[1]])
  p=as.numeric(strsplit(as.character(new_pat$V7[i]), ",")[[1]])
  name$p.value[i] = ks.test(m,p) $ p.value
}

name$p.value.fdr = p.adjust(name$p.value, method = "fdr")

name=name[order(name$p.value, decreasing = FALSE),]

write.table(name, file=paste(Tissue, name_body,"pValue.bed", sep = "_"), quote=F, row.names=F, col.names=F, sep="\t")
#write.bed(name, paste(Tissue, "gencode.vM20.annotation_transcript_30bp_wSNP_10+reads_uniq_pValue.bed", sep = "_"))

