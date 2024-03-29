library("TmCalculator")
library(bigWig)

file_dir="~/Box Sync/BN_IGV/"
setwd("~/Box Sync/Danko_lab_work/F1_8Tissues/Initiation/TSN_ShootingGallery_TSSNotInAlleleHMMBlocks/")

getTSN_2N <- function(seq, d){  
  # generate dinucleotide at position d and d+1 in seq
  seq=s2c(as.character(seq))
  return (c2s(seq[d:(d+1)]))
}

read_read_mat_S <-function (file.plus.bw, file.minus.bw , bed6, step=2, navg = 20, times=1, use.log=FALSE)
{
  bw.plus  <- load.bigWig( file.plus.bw )
  bw.minus <- load.bigWig( file.minus.bw )
  
  hCountMatrix <- bed6.step.bpQuery.bigWig(bw.plus, bw.minus, bed6[,c(1:6)] , step=step, abs.value=TRUE, op = "sum")
  hCountMatrix <- lapply(1:NROW(hCountMatrix), function(i){ if(bed6[i,6]=="-") return(rev(hCountMatrix[[i]])) else return(hCountMatrix[[i]])} );
  if (!use.log){
    hmat <- times * matrix(unlist(hCountMatrix), nrow= NROW(bed6), byrow=TRUE) ;
  } else {
    hmat <- log(times * matrix(unlist(hCountMatrix), nrow= NROW(bed6), byrow=TRUE) + 1) ;
  }
  #avgMat <- t(sapply(1:floor(NROW(hmat)/navg), function(x) {colMeans(hmat[((x-1)*navg+1):min(NROW(hmat),(x*navg)),])}))
  
  unload.bigWig(bw.plus);
  unload.bigWig(bw.minus);
  
  return(hmat);    
}


getInr2N_Dist_Delta_Signal_aroundMaxTSN <- function(df, name, d, WithCAatHighAllele=TRUE, onlyWeakInr=FALSE, onlyCA=FALSE, plot=FALSE){
  if (onlyWeakInr & onlyCA) {
    stop("Invalid 'input' value")
    }
  
  Dist <- NULL
  Delta_Signal <- NULL
  
    cat ("all df ",dim(df), "\n")
    colnames(df)[12:13] = c("HighAlleleSeq", "LowAlleleSeq")
    colnames(df)[9]="Bino_p_value"
    colnames(df)[10] = "strand"
    colnames(df)[11] = "winP"
    
    old_df=df
  for (strand in c("+", "-")){  
    df = old_df[old_df$strand == strand,]
    # keep rows with total allelic reads > 5
    df = df[df$V6+df$V7 >= 5,]
    cat ("df read count >5,  strand= ", strand, "," ,dim(df), "\n")
    
    df$CAatHighTSN =  sapply(df$HighAlleleSeq, getTSN_2N, d=d)
    df$CAatLowTSN =  sapply(df$LowAlleleSeq, getTSN_2N, d=d)
    
    if(WithCAatHighAllele){
      # keep High allele with CA at TSN
      df = df[df$CAatHighTSN=="CA",]
      cat ("df read count >5, High allele with CA, strand= ", strand, "," ,dim(df), "\n")
    }

    # make a matrix of dinucleotide. Sliding windows with overlaps
    h_m <- NULL
    for (i in 1:length(df$HighAlleleSeq)){
      a = s2c(as.character(df$HighAlleleSeq[i]))
      r <- NULL
      for (j in seq(1, 2*d, 1)){
        r = c(r, c2s(a[j:(j+1)]))
      }
      h_m  = rbind(h_m ,r)
    }
    
    l_m <- NULL
    for (i in 1:length(df$LowAlleleSeq)){
      a = s2c(as.character(df$LowAlleleSeq[i]))
      r <- NULL
      for (j in seq(1, 2*d, 1)){
        r = c(r, c2s(a[j:(j+1)]))
      }
      l_m  = rbind(l_m ,r)
    }
    
    dim(h_m)
    dim(l_m)
    
    diNu = c("CA", "CG", "TA", "TG")
    if (onlyWeakInr){
      diNu = c("CG", "TA", "TG")
    }
    if (onlyCA){
      diNu = c("CA")
    }

    
    
    
    # the transcription level near the TSN
    TSS = df[,c(1,2,3,4,5,10)]
    if (strand=="+"){
      TSS[,2] <- df$V2 - d + 1
      TSS[,3] <- df$V3 + d
    }else{
      TSS[,2] <- df$V2 - d 
      TSS[,3] <- df$V3 + d - 1
    }
    
    step=1; times=1; use.log=FALSE
    
    file.bw.plus.pat=paste(file_dir,organ, "_map2ref_1bpbed_map5_CAST_plus.bw", sep="")
    file.bw.minus.pat=paste(file_dir,organ, "_map2ref_1bpbed_map5_CAST_minus.bw", sep="")
    file.bw.plus.mat=paste(file_dir,organ, "_map2ref_1bpbed_map5_B6_plus.bw", sep="")
    file.bw.minus.mat=paste(file_dir,organ, "_map2ref_1bpbed_map5_B6_minus.bw", sep="")
    readCount.pat <- read_read_mat_S (file.bw.plus.pat, file.bw.minus.pat, TSS[,c(1:6)], step, times=times, use.log=use.log)
    readCount.mat <- read_read_mat_S (file.bw.plus.mat, file.bw.minus.mat, TSS[,c(1:6)],  step, times=times, use.log=use.log) 
    readCount.combined <- readCount.pat + readCount.mat
    dim(readCount.pat)
    
    # generate x axis of the plot, distance to the TSN
    x_m  <- NULL
    for (i in 1:dim(readCount.pat)[1]){
      x_m <- rbind( x_m, seq(-d+1,d,1))
    }
    dim(x_m)
    t=1
    # delta_reads =  log2((readCount.HIGH+t)/(readCount.LOW+t))
    delta_reads = log2((readCount.pat+t)/(readCount.mat+t))
    delta_reads[which(df$winP=="M"),] = log2((readCount.mat+t )/( readCount.pat+t))[which(df$winP=="M"),]
    #delta_reads = ((readCount.pat+t)-(readCount.mat+t))
    #delta_reads[which(df$winP=="M"),] = ((readCount.mat+t )-( readCount.pat+t))[which(df$winP=="M"),]
    dim(delta_reads)
    #View(delta_reads)
    #dim(df)
    
    
    # matrix indicates if the dinucleotide is "CA", "CG", "TA", or "TG"
    m = matrix(h_m %in% diNu & readCount.combined >=5, ncol = 100)
    dim(m)
    
    
    if (strand=="+"){
      # plot(x_m[m], delta_reads[m], col=rgb(red=0.2, green=0.2, blue=0.2, alpha=0.15), pch=19, 
      #      #ylim=c(-100,100),
      #      main=name,
      #      ylab = "log2(High Allele + 1 / Low Allele + 1)",
      #      xlab = "Distance to TSN with CA at High Allele")
      Dist <- c(Dist,x_m[m] )
      Delta_Signal <- c(Delta_Signal, delta_reads[m])
      
    }else{
      # points(x_m[m], delta_reads[m], col=rgb(red=0.2, green=0.2, blue=0.2, alpha=0.15), pch=19)
      Dist <- c(Dist, x_m[m])
      Delta_Signal <- c(Delta_Signal, delta_reads[m])
    }
  }
  return (data.frame(dist=Dist, Delta_Signal= Delta_Signal))
}


getInr2N_Delta_Signal_atMaxTSN_HighLow <- function(df, name, d, allele1="CA", allele2="TA"){
  Delta_Signal<-NULL
  cat ("all df ",dim(df), "\n")
  colnames(df)[12:13] = c("HighAlleleSeq", "LowAlleleSeq")
  colnames(df)[9]="Bino_p_value"
  colnames(df)[10] = "strand"
  colnames(df)[11] = "winP"
  
  old_df=df
  for (strand in c("+", "-")){  
    df = old_df[old_df$strand == strand,]
    # keep rows with total allelic reads > 5
    df = df[df$V6+df$V7 >= 5,]
    cat ("df read count >5,  strand= ", strand, "," ,dim(df), "\n")
    
    df$InrAtHighTSN =  sapply(df$HighAlleleSeq, getTSN_2N, d=d)
    df$InrAtLowTSN =  sapply(df$LowAlleleSeq, getTSN_2N, d=d)
    
    df$HighAllele1 = df$InrAtHighTSN == allele1
    df$HighAllele2 = df$InrAtHighTSN == allele2
    df$LowAllele1 = df$InrAtLowTSN == allele1
    df$LowAllele2 = df$InrAtLowTSN == allele2
    
    
    
    # keep those with both allele1 and allele 2
    df = df[(df$HighAllele1|df$LowAllele1)&(df$HighAllele2|df$LowAllele2),]
    cat ("df read count >5, with both allele1 and allele 2, strand= ", strand, "," ,dim(df), "\n")
    if (dim(df)[1]>0) {
    # the transcription level near the TSN
    TSS = df[,c(1,2,3,4,5,10)]
    
    step=1; times=1; use.log=FALSE
    
    file.bw.plus.pat=paste(file_dir,organ, "_map2ref_1bpbed_map5_CAST_plus.bw", sep="")
    file.bw.minus.pat=paste(file_dir,organ, "_map2ref_1bpbed_map5_CAST_minus.bw", sep="")
    file.bw.plus.mat=paste(file_dir,organ, "_map2ref_1bpbed_map5_B6_plus.bw", sep="")
    file.bw.minus.mat=paste(file_dir,organ, "_map2ref_1bpbed_map5_B6_minus.bw", sep="")
    readCount.pat <- read_read_mat_S (file.bw.plus.pat, file.bw.minus.pat, TSS[,c(1:6)], step, times=times, use.log=use.log)
    readCount.mat <- read_read_mat_S (file.bw.plus.mat, file.bw.minus.mat, TSS[,c(1:6)],  step, times=times, use.log=use.log) 
    readCount.combined <- readCount.pat + readCount.mat
    dim(readCount.pat)
    t=1
    df$delta_reads_high_low = ((readCount.pat+t)/(readCount.mat+t))
    df$delta_reads_high_low[which(df$winP=="M"),] = ((readCount.mat+t )/( readCount.pat+t))[which(df$winP=="M"),]
    
    df$delta_reads_allele1_2 = df$delta_reads_high_low
    df$delta_reads_allele1_2[(df$LowAllele1),] = 1/df$delta_reads_high_low[(df$LowAllele1),]
    Delta_Signal <- c(Delta_Signal, log2(df$delta_reads_allele1_2))
    cat("Delta_Signal legnth", length(Delta_Signal), "\n")
    }
    
  }
  return (data.frame(Delta_Signal= Delta_Signal))
}

getInr2N_Delta_Signal_atMaxTSN_matpat <- function(df, name, d, allele1="CA", allele2="TA"){
  Delta_Signal<-NULL
  cat ("all df ",dim(df), "\n")
  colnames(df)[11:12] = c("matAlleleSeq", "patAlleleSeq")
  df$winP = unlist(strsplit(as.character(df$V4), split=","))[seq(1,(length(df$V4)*3)-2,3)]
  colnames(df)[9]="Bino_p_value"
  colnames(df)[10] = "strand"

  old_df=df
  for (strand in c("+", "-")){  
    df = old_df[old_df$strand == strand,]
    # keep rows with total allelic reads > 5
    df = df[df$V6+df$V7 >= 5,]
    cat ("df read count >5,  strand= ", strand, "," ,dim(df), "\n")
    
    df$InrAtMatTSN =  sapply(df$matAlleleSeq, getTSN_2N, d=d)
    df$InrAtPatTSN =  sapply(df$patAlleleSeq, getTSN_2N, d=d)
    
    df$matAllele1 = df$InrAtMatTSN == allele1
    df$matAllele2 = df$InrAtMatTSN == allele2
    df$patAllele1 = df$InrAtPatTSN == allele1
    df$patAllele2 = df$InrAtPatTSN == allele2
    
    
    
    # keep those with both allele1 and allele 2
    df = df[(df$matAllele1|df$patAllele1)&(df$matAllele2|df$patAllele2),]
    cat ("df read count >5, with both allele1 and allele 2, strand= ", strand, "," ,dim(df), "\n")
    if (dim(df)[1]>0) {
      # the transcription level near the TSN
      TSS = df[,c(1,2,3,4,5,10)]
      
      step=1; times=1; use.log=FALSE
      
      file.bw.plus.pat=paste(file_dir,organ, "_map2ref_1bpbed_map5_CAST_plus.bw", sep="")
      file.bw.minus.pat=paste(file_dir,organ, "_map2ref_1bpbed_map5_CAST_minus.bw", sep="")
      file.bw.plus.mat=paste(file_dir,organ, "_map2ref_1bpbed_map5_B6_plus.bw", sep="")
      file.bw.minus.mat=paste(file_dir,organ, "_map2ref_1bpbed_map5_B6_minus.bw", sep="")
      readCount.pat <- read_read_mat_S (file.bw.plus.pat, file.bw.minus.pat, TSS[,c(1:6)], step, times=times, use.log=use.log)
      readCount.mat <- read_read_mat_S (file.bw.plus.mat, file.bw.minus.mat, TSS[,c(1:6)],  step, times=times, use.log=use.log) 
      readCount.combined <- readCount.pat + readCount.mat
      df$readCount.mat = readCount.mat
      df$readCount.pat = readCount.pat

      dim(readCount.pat)
      t=1
      df$delta_reads_mat_pat = ((readCount.mat+t)/(readCount.pat+t))
      #df$delta_reads_high_low[which(df$winP=="M"),] = ((readCount.mat+t )/( readCount.pat+t))[which(df$winP=="M"),]
      
      df$delta_reads_allele1_2 = df$delta_reads_mat_pat 
      df$delta_reads_allele1_2[(df$patAllele1),] = 1/df$delta_reads_mat_pat[(df$patAllele1),]
      Delta_Signal <- c(Delta_Signal, log2(df$delta_reads_allele1_2))
      cat("Delta_Signal legnth", length(Delta_Signal), "\n")
    }
    
  }
  return (data.frame(Delta_Signal= Delta_Signal))
}


# Look at the magnitude of difference between strings in different initiator combinations: 
# CA - vs - TA - vs - TG - vs - CG.
# EXclude TSS within AlleleHMM blocks
setwd("~/Box Sync/Danko_lab_work/F1_8Tissues/Initiation/TSN_ShootingGallery_TSSNotInAlleleHMMBlocks/")
d=50
x_list <- NULL
asTSS=""
SNP_orBackground = "_SNP_"
diNu = c("CA", "TA", "TG", "CG")
#allele1="CA"
#allele2="CG"
#organ="BN"
asTSS_name = asTSS
for (i in 1:3){
  for (j in (i+1):4){
    for (organ in c("BN", "LV")){
      name=paste(organ, asTSS_name, SNP_orBackground, sep = "")
      #df=read.table(file = "BN_allReads_TSS_maxTSNs_SNP_TSSNotInAlleleHMMBlocks_binomtest_+-50_mat_patSeq.bed")
      #df=read.table(paste(organ, "_allReads_TSS_maxTSNs",SNP_orBackground,"TSSNotInAlleleHMMBlocks_binomtest_+-",d,"_High_LowAlleleSeq",asTSS, ".bed", sep=""))
      df=read.table(paste(organ, "_allReads_TSS_maxTSNs",SNP_orBackground,"TSSNotInAlleleHMMBlocks_binomtest_+-",d,"_mat_patSeq",asTSS, ".bed", sep=""))
      cat (paste(organ, "_allReads_TSS_maxTSNs",SNP_orBackground,"TSSNotInAlleleHMMBlocks_binomtest_+-",d,"_mat_patSeq",asTSS, ".bed", sep=""))
      cat ("\n")
      
      cat (i, diNu[i], "\t", j, diNu[j], "\n")
      allele1 = diNu[i]
      allele2 = diNu[j]
      
      temp=getInr2N_Delta_Signal_atMaxTSN_matpat(df, name, d, allele1=allele1, allele2=allele2)
      
      x_list[[paste(name,paste(allele1,allele2, sep="/"), sep="")]] = temp$Delta_Signal
    }
  }
}
str(x_list)





####
# include TSS within AlleleHMM blocks
setwd("~/Box Sync/Danko_lab_work/F1_8Tissues/Initiation/TSN_ShootingGallery_NoAlleleHMMFilter/")
d=50
x_list <- NULL
asTSS=""
SNP_orBackground = "_SNP_"
diNu = c("CA", "TA", "TG", "CG")
#allele1="CA"
#allele2="CG"
#organ="BN"
asTSS_name = asTSS
for (i in 1:3){
  for (j in (i+1):4){
    for (organ in c("BN", "LV")){
      name=paste(organ, asTSS_name, SNP_orBackground, sep = "")
      df=read.table(paste(organ, "_allReads_TSS_maxTSNs",SNP_orBackground,"TSSNoAlleleHMMFilter_binomtest_+-",d,"_mat_patSeq",asTSS, ".bed", sep=""))
      cat (paste(organ, "_allReads_TSS_maxTSNs",SNP_orBackground,"TSSNoAlleleHMMFilter_binomtest_+-",d,"_mat_patSeq",asTSS, ".bed", sep=""))
      cat ("\n")    

      cat (i, diNu[i], "\t", j, diNu[j], "\n")
      allele1 = diNu[i]
      allele2 = diNu[j]
      
      temp=getInr2N_Delta_Signal_atMaxTSN_matpat(df, name, d, allele1=allele1, allele2=allele2)
      
      x_list[[paste(name,paste(allele1,allele2, sep="/"), sep="")]] = temp$Delta_Signal
    }
  }
}
str(x_list)
new_x_list <- NULL
c("CA/TA","CA/TG","CA/CG", "TA/TG")
new_x_list[["CA/TA"]] = c(x_list[["BN_SNP_CA/TA"]],x_list[["LV_SNP_CA/TA"]])
new_x_list[["CA/TG"]] = c(x_list[["BN_SNP_CA/TG"]],x_list[["LV_SNP_CA/TG"]])
new_x_list[["CA/CG"]] = c(x_list[["BN_SNP_CA/CG"]],x_list[["LV_SNP_CA/CG"]])
new_x_list[["TA/TG"]] = c(x_list[["BN_SNP_TA/TG"]],x_list[["LV_SNP_TA/TG"]])
#new_x_list[["TG/CG"]] = c(x_list[["BN_SNP_TG/CG"]],x_list[["LV_SNP_TG/CG"]])

str(new_x_list)
pdf("BNandLV_CA_TA_TG_CG_NoAlleleHMMFilter.pdf", width = 5, height = 5, useDingbats=FALSE)
par(mar=c(5.1, 4.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
library(vioplot)
#par(mgp=c(3,1,0))
vioplot(new_x_list, las=2, 
        #ylim=c(-6,6),
        ylab = "log2(Allele 1 + 1 / Allele 2 + 1)",
        col=c("purple"),
        frame.plot=F
)
#stripchart(new_x_list, vertical = TRUE, add=T, method = "jitter", pch=19, col="black", jitter = 0.1, offset = 0, cex=1, las=2)

abline(h=0, lty=2)
dev.off()


wilcox.test(new_x_list[["CA/TA"]], new_x_list[["CA/TG"]], paired = FALSE)
wilcox.test(new_x_list[["CA/TA"]], new_x_list[["CA/TG"]], paired = FALSE)


wilcox.test(new_x_list[["CA/TA"]], new_x_list[["CA/TG"]], paired = FALSE)
wilcox.test(new_x_list[["CA/TA"]], new_x_list[["CA/CG"]], paired = FALSE)
wilcox.test(new_x_list[["CA/TA"]], new_x_list[["TA/TG"]], paired = FALSE)
wilcox.test(new_x_list[["CA/TG"]], new_x_list[["CA/CG"]], paired = FALSE)
wilcox.test(new_x_list[["CA/TG"]], new_x_list[["TA/TG"]], paired = FALSE)
wilcox.test(new_x_list[["CA/CG"]], new_x_list[["TA/TG"]], paired = FALSE)


wilcox.test(new_x_list[["TA/TG"]])
wilcox.test(new_x_list[["TA/TG"]], rep(0,106), paired = FALSE)



pvalue<-NULL
paires = c("CA/TA","CA/TG","CA/CG", "TA/TG")
name<- NULL
p_value <- NULL
for (i in 1:3){
  for (j in (i+1):4){
    name = c(name, paste(paires[i], paires[j], sep="_"))
    p_value = c(p_value, wilcox.test(new_x_list[[paires[i]]], new_x_list[[paires[j]]], paired = FALSE)$p.value)
  }
}

p.adjust(p_value, method = "fdr")

legend("topleft", legend=c("BN", "LV"),
       #title = "SNPs",
       fill = c("purple", "gray"), 
       bty = "n")


###

###
library("vioplot")
# exclude TSS in AlleleHMM blocks
setwd("~/Box Sync/Danko_lab_work/F1_8Tissues/Initiation/TSN_ShootingGallery_TSSNotInAlleleHMMBlocks/")

d=50
v_list <- NULL
asTSS_list=c("", "_AsTSSSingleBaseDriven", "_AsTSSMultipleBaseDriven")
asTSS_name_list=c("","_Single", "_Multiple")
uper_bound = 20
lower_bound = -20
for (organ in c("BN", "LV")){
  for (a in 1:3){
    asTSS=asTSS_list[a]  ; asTSS
    asTSS_name = asTSS_name_list[a]
    for (SNP_orBackground in c("_NoSNP_","_SNP_")){
      name=paste(organ, asTSS_name, SNP_orBackground, sep = "")
      df=read.table(paste(organ, "_allReads_TSS_maxTSNs",SNP_orBackground,"TSSNotInAlleleHMMBlocks_binomtest_+-",d,"_High_LowAlleleSeq",asTSS, ".bed", sep=""))
      cat (paste(organ, "_allReads_TSS_maxTSNs",SNP_orBackground,"TSSNotInAlleleHMMBlocks_binomtest_+-",d,"_High_LowAlleleSeq",asTSS, ".bed", sep=""))
      cat ("\n")
      temp=getInr2N_Dist_Delta_Signal_aroundMaxTSN(df, name, d)
      temp = temp[temp$dist!=0,]
      v_list[[paste(name, lower_bound, "_",uper_bound, sep="")]] = temp$Delta_Signal[temp$dist >lower_bound & temp$dist<uper_bound ]

      temp=getInr2N_Dist_Delta_Signal_aroundMaxTSN(df, name, d, onlyWeakInr = TRUE)
      temp = temp[temp$dist!=0,]
      v_list[[paste(name, lower_bound,"_",uper_bound,"_OnlyWeakInr", sep="")]] = temp$Delta_Signal[temp$dist >lower_bound & temp$dist<uper_bound ]
      
      temp=getInr2N_Dist_Delta_Signal_aroundMaxTSN(df, name, d, onlyCA = TRUE)
      temp = temp[temp$dist!=0,]
      v_list[[paste(name, lower_bound, "_",uper_bound,"_OnlyCA", sep="")]] = temp$Delta_Signal[temp$dist >lower_bound & temp$dist<uper_bound ]

      str(v_list)
    }
  }
}


SNP_orBackground=c("_SNP_", "_NoSNP_")
SNP_orBackground_name =c("_SNP_","_NoSNP_")
v_list_withNewName <- NULL
for (organ in c("BN", "LV")){
  for (a in 1:3){
    asTSS=asTSS_list[a]
    asTSS_name = asTSS_name_list[a]
    for (k in 1:2){
      name=paste(organ, asTSS_name, sep = "")
      v_list_withNewName[[paste(name,SNP_orBackground_name[k],"AllInr", sep="")]] = v_list[[paste(name,SNP_orBackground[k], lower_bound, "_",uper_bound, sep="")]] 
      v_list_withNewName[[paste(name,SNP_orBackground_name[k],"OnlyWeakInr", sep="")]]=v_list[[paste(name,SNP_orBackground[k], lower_bound,"_",uper_bound,"_OnlyWeakInr", sep="")]]
      v_list_withNewName[[paste(name,SNP_orBackground_name[k],"OnlyCA", sep="")]] = v_list[[paste(name,SNP_orBackground[k], lower_bound, "_",uper_bound,"_OnlyCA", sep="")]] 
    }
  }
}


SNP_orBackground=c("_SNP_", "_NoSNP_")
SNP_orBackground_name =c("_SNP_","_NoSNP_")
Inr_name=c("AllInr","OnlyWeakInr","OnlyCA")
Inr=c("","_OnlyWeakInr","_OnlyCA")
v_list_withNewName <- NULL
w_test <- NULL
for (i in 1:3){
  for (a in 1:3){
    for (organ in c("BN", "LV")){
      asTSS=asTSS_list[a]
      asTSS_name = asTSS_name_list[a]
      for (k in 1:2){
        name=paste(organ, asTSS_name, sep = "")
        v_list_withNewName[[paste(name,SNP_orBackground_name[k],Inr_name[i], sep="")]] = -1 * (v_list[[paste(name,SNP_orBackground[k], lower_bound, "_",uper_bound,Inr[i], sep="")]]) 
        #v_list_withNewName[[paste(name,SNP_orBackground_name[k],"AllInr", sep="")]] = v_list[[paste(name,SNP_orBackground[k], lower_bound, "_",uper_bound, sep="")]] 
        #v_list_withNewName[[paste(name,SNP_orBackground_name[k],"OnlyWeakInr", sep="")]]=v_list[[paste(name,SNP_orBackground[k], lower_bound,"_",uper_bound,"_OnlyWeakInr", sep="")]]
        #v_list_withNewName[[paste(name,SNP_orBackground_name[k],"OnlyCA", sep="")]] = v_list[[paste(name,SNP_orBackground[k], lower_bound, "_",uper_bound,"_OnlyCA", sep="")]] 
      }
     w_test = c(w_test, wilcox.test(v_list_withNewName[[paste(name,SNP_orBackground_name[1],Inr_name[i], sep="")]], v_list_withNewName[[paste(name,SNP_orBackground_name[2],Inr_name[i], sep="")]])$p.value)
    }
  }
}
p.adjust(w_test, method = "fdr")<=0.05
str(v_list_withNewName)
# for sup figure
pdf("ShootingGallery_Inr_comparison.pdf", width = 16, height = 8, useDingbats=FALSE)
par(mar=c(13.1, 4.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
#par(mgp=c(3,1,0))
vioplot(v_list_withNewName, las=2, 
        #ylim=c(-6,6),
        ylab = "log2(Low Allele + 1 / High Allele + 1)",
        #col=c("black", "blue", "orange"), 
        col=c("purple", "gray"),
        #border=c(rep("purple",3) ,rep("gray",3)),
        frame.plot=F)


abline(h=0, lty=2)
legend("topleft", legend=c("maxTSNs with SNPs", "maxTSNs without SNP" ),
       #title = "SNPs",
       fill = c("purple", "gray"), 
       bty = "n")
if(0){
abline(h=0, col="gray")
legend("topleft", legend=c("CA, CG, TA, TG", "CG, TA, TG", "CA only" ),
       title = "Inr",
       fill = c("black", "blue", "orange"),
       bty = "n")
}
dev.off()



# use part of the violin for main figures
new_v_list <- NULL
asTSS_name=""
for (organ in c("BN", "LV")){
  for (SNP_orBackground in c("_SNP_", "_NoSNP_")){
    name=paste(organ, asTSS_name, SNP_orBackground, sep = "")
    new_v_list[[paste(name, sep="")]]=-1*(v_list[[paste(name, lower_bound, "_",uper_bound, sep="")]])
    }
}
str(new_v_list)

pdf("ShootingGallery.pdf", width = 5, height = 5, useDingbats=FALSE)
par(mar=c(5.1, 4.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1

#par(mgp=c(3,1,0))
vioplot(new_v_list, las=1, 
        #ylim=c(-6,6),
        ylab="Log2(Non-CA allele / CA allele)",
        #ylab = "log2(Low Allele + 1 / High Allele + 1)",
        col=c("purple", "gray"),
        frame.plot=F
        )

abline(h=0, lty=2)
legend("topleft", legend=c("maxTSNs with SNPs", "maxTSNs without SNP"),
       #title = "SNPs",
       fill = c("purple", "gray"), 
       bty = "n")

dev.off()


# Wilcoxon Rank Sum and Signed Rank Tests

wilcox.test(new_v_list$BN_SNP_, new_v_list$BN_NoSNP_)
wilcox.test(new_v_list$LV_SNP_, new_v_list$LV_NoSNP_)

wilcox.test(v_list_withNewName$BN_Single_SNP_AllInr, v_list_withNewName$BN_Single_NoSNP_AllInr)





# scatter plots of BN
d=50
asTSS = ""
asTSS_name=asTSS
organ="BN"
SNP_orBackground ="_NoSNP_"

name=paste(organ, asTSS_name, SNP_orBackground, sep = "")
df=read.table(paste(organ, "_allReads_TSS_maxTSNs",SNP_orBackground,"TSSNotInAlleleHMMBlocks_binomtest_+-",d,"_High_LowAlleleSeq",asTSS, ".bed", sep=""))
temp=getInr2N_Dist_Delta_Signal_aroundMaxTSN(df, name, d)

SNP_orBackground ="_SNP_"

name=paste(organ, asTSS_name, SNP_orBackground, sep = "")
df=read.table(paste(organ, "_allReads_TSS_maxTSNs",SNP_orBackground,"TSSNotInAlleleHMMBlocks_binomtest_+-",d,"_High_LowAlleleSeq",asTSS, ".bed", sep=""))
temp2=getInr2N_Dist_Delta_Signal_aroundMaxTSN(df, name, d)
out=temp2


pdf("ShootingGallery_scatter_BN-2.pdf", width = 8, height = 8, useDingbats=FALSE)
par(mfcol=c(3,2))

#out=temp2
name=""
for (out in list(temp, temp2)){
plot(out$dist, out$Delta_Signal, col=rgb(red=0.2, green=0.2, blue=0.2, alpha=0.15), pch=19,
     #ylim=c(-100,100),
     main=name,
     ylab = "log2(High Allele + 1 / Low Allele + 1)",
     xlab = "Distance to maxTSN with CA at High Allele",
     las=1,
     #ylim=c(-6,6),
     frame.plot=F)
  abline(v=20, lty=2, col="blue")
  abline(v=-20, lty=2, col="blue")
  abline(h=0)

out = out[out$dist!=0,]

signal.lo <-  loess(Delta_Signal ~ dist, out)
x=seq(-d+1,d,1)
p=predict(signal.lo, data.frame(dist = x), se = TRUE)
plot(out, col=rgb(red=0.2, green=0.2, blue=0.2, alpha=0.15), pch=19,
     main=name,
     #ylim=c(-5,5),
     ylab = "log2(High Allele + 1 / Low Allele + 1)",
     xlab = "Distance to maxTSN with CA at High Allele", 
     #ylim=c(-4,4)
     las=1,
     frame.plot=F)

abline(v=20, lty=2, col="blue")
abline(v=-20, lty=2, col="blue")
abline(h=0)
lines(x, p$fit, col="red" )
#abline(v=0)
plot(x, p$fit, col="red",  ylab = "Predicted log2(High Allele + 1 / Low Allele + 1)",
     main=name,
     ylim=c(-0.5,0.5),
     xlab = "Distance to TSN with CA at High Allele",
     las=1, type="l",
     frame.plot=F, lwd=3)
abline(h=0)
#abline(v=0)
abline(v=20, lty=2, col="blue")
abline(v=-20, lty=2, col="blue")

}


dev.off()


# scatter plots of LV
d=50
asTSS = ""
asTSS_name=asTSS
organ="LV"
SNP_orBackground ="_NoSNP_"

name=paste(organ, asTSS_name, SNP_orBackground, sep = "")
df=read.table(paste(organ, "_allReads_TSS_maxTSNs",SNP_orBackground,"TSSNotInAlleleHMMBlocks_binomtest_+-",d,"_High_LowAlleleSeq",asTSS, ".bed", sep=""))
temp=getInr2N_Dist_Delta_Signal_aroundMaxTSN(df, name, d)

SNP_orBackground ="_SNP_"

name=paste(organ, asTSS_name, SNP_orBackground, sep = "")
df=read.table(paste(organ, "_allReads_TSS_maxTSNs",SNP_orBackground,"TSSNotInAlleleHMMBlocks_binomtest_+-",d,"_High_LowAlleleSeq",asTSS, ".bed", sep=""))
temp2=getInr2N_Dist_Delta_Signal_aroundMaxTSN(df, name, d)
out=temp2


pdf("ShootingGallery_scatter_LV.pdf", width = 8, height = 8, useDingbats=FALSE)
par(mfcol=c(3,2))

#out=temp2
name=""
for (out in list(temp, temp2)){
  plot(out$dist, out$Delta_Signal, col=rgb(red=0.2, green=0.2, blue=0.2, alpha=0.15), pch=19,
       #ylim=c(-100,100),
       main=name,
       ylab = "log2(High Allele + 1 / Low Allele + 1)",
       xlab = "Distance to maxTSN with CA at High Allele",
       las=1,
       #ylim=c(-6,6),
       frame.plot=F)
  abline(v=20, lty=2, col="blue")
  abline(v=-20, lty=2, col="blue")
  abline(h=0)
  
  out = out[out$dist!=0,]
  
  signal.lo <-  loess(Delta_Signal ~ dist, out)
  x=seq(-d+1,d,1)
  p=predict(signal.lo, data.frame(dist = x), se = TRUE)
  plot(out, col=rgb(red=0.2, green=0.2, blue=0.2, alpha=0.15), pch=19,
       main=name,
       #ylim=c(-5,5),
       ylab = "log2(High Allele + 1 / Low Allele + 1)",
       xlab = "Distance to maxTSN with CA at High Allele", 
       #ylim=c(-4,4)
       las=1,
       frame.plot=F)
  
  abline(v=20, lty=2, col="blue")
  abline(v=-20, lty=2, col="blue")
  abline(h=0)
  lines(x, p$fit, col="red" )
  #abline(v=0)
  plot(x, p$fit, col="red",  ylab = "Predicted log2(High Allele + 1 / Low Allele + 1)",
       main=name,
       ylim=c(-0.5,0.5),
       xlab = "Distance to TSN with CA at High Allele",
       las=1, type="l",
       frame.plot=F, lwd=3)
  abline(h=0)
  #abline(v=0)
  abline(v=20, lty=2, col="blue")
  abline(v=-20, lty=2, col="blue")
  
}


dev.off()


d=50

S_list <- NULL
asTSS=""
asTSS_name=""
for (organ in c("BN", "LV")){
  for (SNP_orBackground in c("_NoSNP_","_SNP_")){
    
    name=paste(organ, asTSS_name, SNP_orBackground, sep = "")
    # exclude TSS overlap with AlleleHMM blocks
    df=read.table(paste(organ, "_allReads_TSS_maxTSNs",SNP_orBackground,"TSSNotInAlleleHMMBlocks_binomtest_+-",d,"_High_LowAlleleSeq",asTSS, ".bed", sep=""))
    # include TSS overlap with AlleleHMM blocks
    #df=read.table(paste(organ, "_allReads_TSS_maxTSNs",SNP_orBackground,"binomtest_+-",d,"_High_LowAlleleSeq",asTSS, ".bed", sep=""))
    temp=getInr2N_Dist_Delta_Signal_aroundMaxTSN(df, name, d)
    temp = temp[temp$dist!=0,]
    S_list[[name]]=temp
  }
}
str(S_list)
organ="BN"
organ="LV"
b_list <- NULL
step=5
SNP_orBackground=c("_SNP_", "_NoSNP_")
SNP_orBackground_name =c("_SNP_","_NoSNP_")
for (lower_bound in seq(-40,30,step)){
  for (s in 1:2){
      uper_bound = lower_bound + step
      name=paste(organ, asTSS_name, SNP_orBackground[s], sep = "")
      new_name = paste(organ, asTSS_name, SNP_orBackground_name[s], sep = "")
      temp=S_list[[name]]
      b_list[[paste(new_name, lower_bound, "_",uper_bound, sep="")]] = temp$Delta_Signal[temp$dist >= lower_bound & temp$dist < uper_bound ]
      
      
      str(b_list)
    }
  }

# combine BN and LV
organ="BN+LV"

b_list <- NULL
step=5
SNP_orBackground=c("_SNP_", "_NoSNP_")
SNP_orBackground_name =c("_SNP_","_NoSNP_")
for (lower_bound in seq(-25,20,step)){
  for (s in 1:2){
    uper_bound = lower_bound + step
    new_name = paste(organ, asTSS_name, SNP_orBackground_name[s], sep = "")
    name1 = paste("BN", asTSS_name, SNP_orBackground[s], sep = "")
    name2 = paste("LV", asTSS_name, SNP_orBackground[s], sep = "")
    temp1=S_list[[name1]]
    temp2=S_list[[name2]]
    b_list[[paste(new_name, lower_bound, "_",uper_bound, sep="")]] = c(-1*(temp1$Delta_Signal[temp1$dist >= lower_bound & temp1$dist < uper_bound ]),
                                                                       -1*(temp2$Delta_Signal[temp2$dist >= lower_bound & temp2$dist < uper_bound ]))
    
    
    str(b_list)
  }
}





pdf("ShootingGallery_boxplot_BN+LV_AlleleHMMoverlappingTSSRemoved.pdf", width = 10, height = 6, useDingbats=FALSE)
par(mar=c(10.1, 4.1, 2.1, 2.1))
at.x <- NULL # set here the X-axis positions
x=1
for (i in 1:(length(b_list)/2)){
  at.x <- c(at.x,x)
  x = x+1
  at.x <- c(at.x,x)
  x = x+1.5
}

boxplot(b_list, las=2, 
        col=c("purple","gray"),
        ylab = "log2(Low Allele + 1 / High Allele + 1)",
        frame.plot=F,
        outline=FALSE,
        #space = rep(c(1,2),(length(b_list)/2))
        at=at.x
        )
abline(h=0, lty=2)

legend("topleft", legend=c("maxTSNs with SNPs", "maxTSNs without SNP"),
       #title = "SNPs",
       fill = c("purple", "gray"), 
       bty = "n")
dev.off()
# Wilcoxon Rank Sum and Signed Rank Tests

wilcox.test(new_v_list$BN_SNP_, new_v_list$BN_NoSNP_)
wilcox.test(new_v_list$LV_SNP_, new_v_list$LV_NoSNP_)

w_p_value <- NULL
for (lower_bound in seq(-25,20,step)){
  uper_bound = lower_bound + step
  name1 = paste(organ, asTSS_name, "_SNP_",lower_bound, "_",uper_bound, sep = "")
  name2 = paste(organ, asTSS_name, "_NoSNP_",lower_bound, "_",uper_bound, sep = "")
  w_p_value= c(w_p_value, wilcox.test(b_list[[name1]], b_list[[name2]])$p.value)
#str(w)
}

p.adjust(w_p_value, method = "fdr")



# as a supp figs
# ShootingGallery_boxplot_BN+LV_WithoutAlleleHMMFilter
# include TSS overlapping with AlleleHMM blocks
setwd("~/Box Sync/Danko_lab_work/F1_8Tissues/Initiation/TSN_ShootingGallery_NoAlleleHMMFilter/")

d=50

S_list <- NULL
asTSS=""
asTSS_name=""
for (organ in c("BN", "LV")){
  for (SNP_orBackground in c("_NoSNP_","_SNP_")){
    
    name=paste(organ, asTSS_name, SNP_orBackground, sep = "")
    # exclude TSS overlap with AlleleHMM blocks
    df=read.table(paste(organ, "_allReads_TSS_maxTSNs",SNP_orBackground,"TSSNoAlleleHMMFilter_binomtest_+-",d,"_High_LowAlleleSeq",asTSS, ".bed", sep=""))
    # include TSS overlap with AlleleHMM blocks
    #df=read.table(paste(organ, "_allReads_TSS_maxTSNs",SNP_orBackground,"binomtest_+-",d,"_High_LowAlleleSeq",asTSS, ".bed", sep=""))
    temp=getInr2N_Dist_Delta_Signal_aroundMaxTSN(df, name, d)
    temp = temp[temp$dist!=0,]
    S_list[[name]]=temp
  }
}

str(S_list)

# combine BN and LV
organ="BN+LV"

b_list <- NULL
step=5
SNP_orBackground=c("_SNP_", "_NoSNP_")
SNP_orBackground_name =c("_SNP_","_NoSNP_")
for (lower_bound in seq(-25,20,step)){
  for (s in 1:2){
    uper_bound = lower_bound + step
    new_name = paste(organ, asTSS_name, SNP_orBackground_name[s], sep = "")
    name1 = paste("BN", asTSS_name, SNP_orBackground[s], sep = "")
    name2 = paste("LV", asTSS_name, SNP_orBackground[s], sep = "")
    temp1=S_list[[name1]]
    temp2=S_list[[name2]]
    b_list[[paste(new_name, lower_bound, "_",uper_bound, sep="")]] = c(temp1$Delta_Signal[temp1$dist >= lower_bound & temp1$dist < uper_bound ],
                                                                       temp2$Delta_Signal[temp2$dist >= lower_bound & temp2$dist < uper_bound ])
    
    
    str(b_list)
  }
}


pdf("ShootingGallery_boxplot_BN+LV_WithoutAlleleHMMFilter.pdf", width = 10, height = 6, useDingbats=FALSE)
par(mar=c(10.1, 4.1, 2.1, 2.1))
at.x <- NULL # set here the X-axis positions
x=1
for (i in 1:(length(b_list)/2)){
  at.x <- c(at.x,x)
  x = x+1
  at.x <- c(at.x,x)
  x = x+1.5
}

boxplot(b_list, las=2, 
        col=c("purple","gray"),
        ylab = "log2(High Allele + 1 / Low Allele + 1)",
        frame.plot=F,
        outline=FALSE,
        #space = rep(c(1,2),(length(b_list)/2))
        at=at.x
)
abline(h=0, lty=2)

legend("topright", legend=c("maxTSNs with SNPs", "maxTSNs without SNP"),
       #title = "SNPs",
       fill = c("purple", "gray"), 
       bty = "n")
dev.off
# Wilcoxon Rank Sum and Signed Rank Tests

w_p_value <- NULL
for (lower_bound in seq(-25,20,step)){
  uper_bound = lower_bound + step
  name1 = paste(organ, asTSS_name, "_SNP_",lower_bound, "_",uper_bound, sep = "")
  name2 = paste(organ, asTSS_name, "_NoSNP_",lower_bound, "_",uper_bound, sep = "")
  w_p_value= c(w_p_value, wilcox.test(b_list[[name1]], b_list[[name2]])$p.value)
  #str(w)
}

p.adjust(w_p_value, method = "fdr")



