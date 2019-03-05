#R --vanilla --slave --args $(pwd) PREFIX < sen_spec_matrix.R

#arguments here
args=(commandArgs(TRUE))
setwd(args[1])
PREFIX=args[2]

#setwd("~/Box Sync/Danko_lab_work/F1_8Tissues/sen_spec_matrix")

library("RColorBrewer")
red_cols = brewer.pal(n = 7, name = "Reds")
blue_cols = brewer.pal(n = 7, name = "Blues")

AlleleHMM_REGION_sensitivity_Specificity <- function (input_fp) {
  df=read.table(input_fp)
  df=df[df$V12+df$V13 >=10,]  # reads count within gene >=10
  colnames(df)[7]="GeneAS"
  colnames(df)[9]="AlleleHMMAS"
  colnames(df)[8]="AlleleDBAS"
  diff=0.05
  geneBiased=df[df$GeneAS!="S" & (df$V10<=0.2 |df$V10 >= 0.8 ),]
  geneSym=df[df$GeneAS=="S"& (df$V10>0.45 & df$V10 < 0.55 ),]
  
  df=rbind.data.frame(geneBiased, geneSym)
  AlleleHMM_result = data.frame(mat=seq(diff,1,diff), S_S=0, S_P=0, S_M=0, P_S=0, P_P=0, P_M=0, M_S=0, M_P=0, M_M=0)
  AS=c("S","P","M")
  for (m in 1:(dim(AlleleHMM_result)[1])){
    sub_df = df[df$V10 <= AlleleHMM_result$mat[m] & df$V10 >= (AlleleHMM_result$mat[m] - diff) ,]
    #V10 = mat/mat+pat
    for (i in seq(3)){
      for (j in seq(3)){
        AlleleHMM_result[m,3*(i-1)+j+1]=sum(sub_df$GeneAS==AS[i] & sub_df$AlleleHMMAS==AS[j])
      }
    }
  }
  AlleleHMM_result$mat = AlleleHMM_result$mat - diff/2
  AlleleHMM_result$sen = with(AlleleHMM_result, (sum(P_P)+sum(M_M))/(sum(P_S)+sum(P_P)+sum(M_M)+sum(M_S)))
  AlleleHMM_result$spec = with(AlleleHMM_result, sum(S_S)/(sum(S_S)+sum(S_P)+sum(S_M)))
  AlleleHMM_result$prec = with(AlleleHMM_result, (sum(P_P)+sum(M_M))/(sum(P_P)+sum(M_M)+sum(S_P)+sum(S_M)))
  return (list(AlleleHMM_result$mat, AlleleHMM_result$sen[1], AlleleHMM_result$spec[1], AlleleHMM_result$prec[1]))
}

AlleleDB_REGION_sensitivity_Specificity <- function (input_fp) {
  df=read.table(input_fp)
  df=df[df$V12+df$V13 >=10,]
  colnames(df)[7]="GeneAS"
  colnames(df)[9]="AlleleHMMAS"
  colnames(df)[8]="AlleleDBAS"
  diff=0.05
  geneSym=df[df$GeneAS=="S",]
  geneBiased=df[df$GeneAS!="S" & (df$V10<=0.2 |df$V10 >= 0.8 ),]
  df=rbind.data.frame(geneBiased, geneSym)
  AlleleDB_result = data.frame(mat=seq(diff,1,diff), S_S=0, S_P=0, S_M=0, P_S=0, P_P=0, P_M=0, M_S=0, M_P=0, M_M=0)
  AS=c("S","P","M")
  for (m in 1:(dim(AlleleDB_result)[1])){
    sub_df = df[df$V10 <= AlleleDB_result$mat[m] & df$V10 >= (AlleleDB_result$mat[m] - diff) ,]
    
    for (i in seq(3)){
      for (j in seq(3)){
        AlleleDB_result[m,3*(i-1)+j+1]=sum(sub_df$GeneAS==AS[i] & sub_df$AlleleDBAS==AS[j])
      }
    }
  }
  AlleleDB_result$mat = AlleleDB_result$mat - diff/2
  AlleleDB_result$sen = with(AlleleDB_result, (sum(P_P)+sum(M_M))/(sum(P_S)+sum(P_P)+sum(M_M)+sum(M_S)))
  AlleleDB_result$spec = with(AlleleDB_result, sum(S_S)/(sum(S_S)+sum(S_P)+sum(S_M)))
  AlleleDB_result$prec = with(AlleleDB_result, (sum(P_P)+sum(M_M))/(sum(P_P)+sum(M_M)+sum(S_P)+sum(S_M)))
  return (list(AlleleDB_result$mat, AlleleDB_result$sen[1], AlleleDB_result$spec[1], AlleleDB_result$prec[1]))
}


#PREFIX="BN_PB6_all_R1"
pdf(paste(PREFIX,"_tao_performance_geneSym_and_geneBiased28.pdf",sep=""))
par(mar=c(6.1, 7.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
par(mgp=c(3,1,0))
par(cex.lab=2.2, cex.axis=2.2)
par(bty = 'n') 
i=9
a=AlleleHMM_REGION_sensitivity_Specificity(paste(PREFIX,"_SNP_Allele_Specificity_matrix_t1E-0",i,".bed",sep=""))
d=AlleleDB_REGION_sensitivity_Specificity(paste(PREFIX,"_SNP_Allele_Specificity_matrix_t1E-0",i,".bed",sep=""))
tao=paste("1E-0",i,sep="")
sen=a[[2]]
spec=a[[3]]
prec=a[[4]]
for (i in 8:1){
  a=AlleleHMM_REGION_sensitivity_Specificity(paste(PREFIX,"_SNP_Allele_Specificity_matrix_t1E-0",i,".bed",sep=""))
  tao=c(tao,paste("1E-0",i,sep=""))
  sen=c(sen, a[[2]])
  spec=c(spec, a[[3]])
  prec=c(prec,a[[4]])
}
plot(tao, sen, log="x", type='b', pch=1, ylab="Value", col="red", ylim=c(0,1), xlab="Tuning parameter Tao", cex=1.5,lwd=2, las=1)
lines(tao, spec, type='b', pch=0, lty=1, col="blue", ylim=c(0,1), main=PREFIX, xlab="Maternal reads fraction", cex=1.5,lwd=2, las=1)
lines(tao, prec, type='b', pch=2, lty=1, col="green", ylim=c(0,1), main=PREFIX, xlab="Maternal reads fraction", cex=1.5,lwd=2, las=1)

legend("bottomright", lty =1,col=c("blue", "red", "green"),pch=c(0,1,2),
       legend=c("Specificity", "Sensitivity", "Precision"), cex=1.5,lwd=2, bty="n", title = PREFIX)
dev.off()



AlleleHMM_sensitivity_Specificity <- function (input_fp) {
  df=read.table(input_fp)
  df=df[df$V12+df$V13 >=10,]
  colnames(df)[7]="GeneAS"
  colnames(df)[9]="AlleleHMMAS"
  colnames(df)[8]="AlleleDBAS"
  diff=0.05
  AlleleHMM_result = data.frame(mat=seq(diff,1,diff), S_S=0, S_P=0, S_M=0, P_S=0, P_P=0, P_M=0, M_S=0, M_P=0, M_M=0)
  AS=c("S","P","M")
  for (m in 1:(dim(AlleleHMM_result)[1])){
    sub_df = df[df$V10 <= AlleleHMM_result$mat[m] & df$V10 >= (AlleleHMM_result$mat[m] - diff) ,]
    #V10 = mat/mat+pat
    for (i in seq(3)){
      for (j in seq(3)){
        AlleleHMM_result[m,3*(i-1)+j+1]=sum(sub_df$GeneAS==AS[i] & sub_df$AlleleHMMAS==AS[j])
      }
    }
  }
  AlleleHMM_result$mat = AlleleHMM_result$mat - diff/2
  AlleleHMM_result$sen = with(AlleleHMM_result, (P_P+M_M)/(P_S+P_P+M_M+M_S))
  AlleleHMM_result$spec = with(AlleleHMM_result, S_S/(S_S+S_P+S_M))
  AlleleHMM_result$prec = with(AlleleHMM_result, (P_P+M_M)/(P_P+M_M+S_M+S_P))
  return (list(AlleleHMM_result$mat, AlleleHMM_result$sen, AlleleHMM_result$spec, AlleleHMM_result$prec))
}

AlleleDB_sensitivity_Specificity <- function (input_fp) {
  df=read.table(input_fp)
  df=df[df$V12+df$V13 >=10,]
  colnames(df)[7]="GeneAS"
  colnames(df)[9]="AlleleHMMAS"
  colnames(df)[8]="AlleleDBAS"
  diff=0.05
  AlleleDB_result = data.frame(mat=seq(diff,1,diff), S_S=0, S_P=0, S_M=0, P_S=0, P_P=0, P_M=0, M_S=0, M_P=0, M_M=0)
  AS=c("S","P","M")
  for (m in 1:(dim(AlleleDB_result)[1])){
    sub_df = df[df$V10 <= AlleleDB_result$mat[m] & df$V10 >= (AlleleDB_result$mat[m] - diff) ,]
    
    for (i in seq(3)){
      for (j in seq(3)){
        AlleleDB_result[m,3*(i-1)+j+1]=sum(sub_df$GeneAS==AS[i] & sub_df$AlleleDBAS==AS[j])
      }
    }
  }
  AlleleDB_result$mat = AlleleDB_result$mat - diff/2
  AlleleDB_result$sen = with(AlleleDB_result, (P_P+M_M)/(P_S+P_P+M_M+M_S))
  AlleleDB_result$spec = with(AlleleDB_result, S_S/(S_S+S_P+S_M))
  AlleleDB_result$prec = with(AlleleDB_result, (P_P+M_M)/(P_P+M_M+S_M+S_P))
  return (list(AlleleDB_result$mat, AlleleDB_result$sen, AlleleDB_result$spec, AlleleDB_result$prec))
}


pdf(paste(PREFIX,"_sensitivity_tao_performance.pdf",sep=""))
par(mar=c(6.1, 7.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
par(mgp=c(3,1,0))
par(cex.lab=2.2, cex.axis=2.2)
par(bty = 'n') 
#display.brewer.pal(n = 5, name = "OrRd")
red_cols = brewer.pal(n = 9, name = "OrRd")
type='b'
i=9
a=AlleleHMM_sensitivity_Specificity(paste(PREFIX,"_SNP_Allele_Specificity_matrix_t1E-0",i,".bed",sep=""))
plot(a[[1]], a[[2]], type=type, pch=20, lty=1, ylab="Sensitivity", col=red_cols[i], ylim=c(0,1), xlab="Maternal reads fraction", cex=1.5,lwd=2, las=1)
d=AlleleDB_sensitivity_Specificity(paste(PREFIX,"_SNP_Allele_Specificity_matrix_t1E-0",i,".bed",sep=""))
lines(d[[1]], d[[2]], type=type, pch=20, lty=1, col="blue", cex=1.5,lwd=2)
for (i in 8:1){
  a=AlleleHMM_sensitivity_Specificity(paste(PREFIX,"_SNP_Allele_Specificity_matrix_t1E-0",i,".bed",sep=""))
  lines(a[[1]], a[[2]],type=type, pch=20, lty=(i%%2+1), ylab="Sensitivity", col=red_cols[i], cex=1.5,lwd=2)
}
dev.off()


pdf(paste(PREFIX,"_precision_tao_performance.pdf",sep=""))
par(mar=c(6.1, 7.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
par(mgp=c(3,1,0))
par(cex.lab=2.2, cex.axis=2.2)
par(bty = 'n') 
#display.brewer.pal(n = 9, name = "OrRd")
red_cols = brewer.pal(n = 9, name = "OrRd")
type='b'
i=9
a=AlleleHMM_sensitivity_Specificity(paste(PREFIX,"_SNP_Allele_Specificity_matrix_t1E-0",i,".bed",sep=""))
plot(a[[1]], a[[4]], type=type, pch=20, lty=1, ylab="Precision", col=red_cols[i], ylim=c(0,1), xlab="Maternal reads fraction", cex=1.5,lwd=2, las=1)
d=AlleleDB_sensitivity_Specificity(paste(PREFIX,"_SNP_Allele_Specificity_matrix_t1E-0",i,".bed",sep=""))
lines(d[[1]], d[[4]], type=type, pch=20, lty=1, col="blue", cex=1.5,lwd=2)
for (i in 8:1){
  a=AlleleHMM_sensitivity_Specificity(paste(PREFIX,"_SNP_Allele_Specificity_matrix_t1E-0",i,".bed",sep=""))
  lines(a[[1]], a[[4]],type=type, pch=20, lty=(i%%2+1), col=red_cols[i], cex=1.5,lwd=2)
}
dev.off()
