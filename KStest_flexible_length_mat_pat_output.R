#R --vanilla --slave --args $(pwd) mat_name pat_name output < KStest.R
# only use bed regions with at least 5 VALID mat AND 5 pat VALID reads
# input file format
# chr9    121858088      121858089    19,17       111 +   16,16,16,16,16,18,18,18,18,18,18,19,19,21,21,21,21
# chr9    122805884      122805885    19,25       111 +   17,17,17,17,17,17,17,17,17,17,17,18,18,20,21,21,21,21,21,21,21,21,21,21,21
# chr9    122950964      122950965    20,38       111 -   35,24,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,20,20,20,19,19,19,19,19,19,19,19,16,16,16,16,16
# chr9    123366949      123366950    104,54      111 +   16,20,22,22,22,22,22,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,24,24,24,24,24,24,24,24,24,24,24,24,25,25,25,25,25,25,26,26,27,27,28,29,29,29,30,30,30,31
# chr9    123851874      123851875    19,13       111 -   28,25,23,23,23,22,21,21,20,20,20,20,16
# chr9    123852160      123852161    15,22       111 +   15,15,15,15,18,18,18,18,18,18,19,21,21,21,21,21,21,21,21,21,21,21


#arguments here
args=(commandArgs(TRUE))
setwd(args[1])
mat_name=args[2]
pat_name=args[3]
output=args[4]
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


mat<-read.table(mat_name)
pat<-read.table(pat_name)

for (i in 1:dim(mat)[1]){
  mat$read.count[i] = length(as.numeric(strsplit(as.character(mat$V7[i]), ",")[[1]]))
  pat$read.count[i] = length(as.numeric(strsplit(as.character(pat$V7[i]), ",")[[1]]))
  }

new_mat = mat[(mat$read.count>=5 & pat$read.count>=5),]
new_pat = pat[(mat$read.count>=5 & pat$read.count>=5),]
name = new_mat[1:6]
name$mapped.read.count = paste(new_mat$read.count, new_pat$read.count, sep=",")


for (i in 1:dim(new_mat)[1]){
  m=as.numeric(strsplit(as.character(new_mat$V7[i]), ",")[[1]])
  p=as.numeric(strsplit(as.character(new_pat$V7[i]), ",")[[1]])
  name$p.value[i] = ks.test(m,p) $ p.value
}

name$p.value.fdr = p.adjust(name$p.value, method = "fdr")

name=name[order(name$p.value, decreasing = FALSE),]

write.table(name, file=output, quote=F, row.names=F, col.names=F, sep="\t")

