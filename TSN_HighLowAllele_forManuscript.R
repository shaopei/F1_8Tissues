library("TmCalculator")
library(bigWig)

file_dir="~/Box Sync/BN_IGV/"
setwd("~/Box Sync/Danko_lab_work/F1_8Tissues/Initiation/TSN_ShootingGallery/")

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


getInr2N_Dist_Delta_Signal_aroundMaxTSN <- function(df, name, d, WithCAatHighAllele=TRUE, onlyWeakInr=FALSE, onlyCA=FALSE){
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

organ="BN"; 
#asTSS ="SingleBaseDriven" 
#asTSS="MultipleBaseDriven"
asTSS=""
#SNP_orBackground=""
SNP_orBackground="_SNP"
d=50; name=paste(organ, asTSS, SNP_orBackground, sep = " ")




# use all TSN with at least 5 reads (allelic or not)
#df=read.table(paste(organ, "_allReads_TSN5+_SNP_binomtest_interestingHets_+-",d,"_High_LowAlleleSeq.bed", sep=""))
#df=read.table(paste(organ, "_allReads_TSN5+_SNP_binomtest_+-",d,"_High_LowAlleleSeq.bed", sep=""))

# use maxTSNs only
df=read.table(paste(organ, "_allReads_TSS_maxTSNs",SNP_orBackground,"_TSSNotInAlleleHMMBlocks_binomtest_+-",d,"_High_LowAlleleSeq.bed", sep=""))

# use maxTSNs inside AsTSS, single or multiple base driven
#df=read.table(paste(organ, "_allReads_TSS_maxTSNs_SNP_TSSNotInAlleleHMMBlocks_binomtest_interestingHets_+-",d,"_High_LowAlleleSeq.bed", sep=""))
#df=read.table(paste(organ, "_allReads_TSS_maxTSNs_SNP_TSSNotInAlleleHMMBlocks_binomtest_+-",d,"_High_LowAlleleSeq_AsTSS",asTSS, ".bed", sep=""))

raw_out = getInr2N_Dist_Delta_Signal_aroundMaxTSN(df, name, d)

par(mfrow=c(3,1))

plot(raw_out, col=rgb(red=0.2, green=0.2, blue=0.2, alpha=0.15), pch=19,
     main=name,
     ylab = "log2(High Allele + 1 / Low Allele + 1)",
     xlab = "Distance to TSN with CA at High Allele", 
     #ylim=c(-4,4)
)
abline(h=0)
#out = data.frame(dist=Dist, y= Delta_Signal)
out = raw_out[raw_out $dist!=0,]
#out = out[out$y<0 , ]
signal.lo <-  loess(Delta_Signal ~ dist, out)
x=seq(-d+1,d,1)
p=predict(signal.lo, data.frame(dist = x), se = TRUE)
plot(out, col=rgb(red=0.2, green=0.2, blue=0.2, alpha=0.15), pch=19,
     main=name,
     ylab = "log2(High Allele + 1 / Low Allele + 1)",
     xlab = "Distance to TSN with CA at High Allele", 
     #ylim=c(-4,4)
     )
abline(h=0)
lines(x, p$fit, col="red" )
abline(v=0)
plot(x, p$fit, col="red",  ylab = "Predicted log2(High Allele + 1 / Low Allele + 1)",
     main=name,
     #ylim=c(-0.15, 0.01),
     xlab = "Distance to TSN with CA at High Allele")
abline(h=0)
abline(v=0)

###
library("vioplot")


d=50
v_list <- NULL
#v_name_list <- NULL
asTSS_list=c("_AsTSSSingleBaseDriven", "_AsTSSMultipleBaseDriven", "")
asTSS_name_list=c("_Single", "_Multiple", "")
uper_bound = 20
lower_bound = -20
for (organ in c("BN", "LV")){
  for (a in 1:3){
    asTSS=asTSS_list[a]
    asTSS_name = asTSS_name_list[a]
    for (SNP_orBackground in c("_","_SNP_")){
      name=paste(organ, asTSS_name, SNP_orBackground, sep = "")
      df=read.table(paste(organ, "_allReads_TSS_maxTSNs",SNP_orBackground,"TSSNotInAlleleHMMBlocks_binomtest_+-",d,"_High_LowAlleleSeq",asTSS, ".bed", sep=""))
      temp=getInr2N_Dist_Delta_Signal_aroundMaxTSN(df, name, d)
      temp = temp[temp$dist!=0,]
      v_list[[paste(name, "-20_",uper_bound, sep="")]] = temp$Delta_Signal[temp$dist >-20 & temp$dist<uper_bound ]
      #v_name_list <- c(v_name_list, paste(name, "-20_0", sep="") )
      # v_list[[paste(name, "0_20", sep="")]] = temp$Delta_Signal[temp$dist >0 & temp$dist<20 ]
      # v_name_list <- c(v_name_list, paste(name, "0_20", sep="") )
      
      temp=getInr2N_Dist_Delta_Signal_aroundMaxTSN(df, name, d, onlyWeakInr = TRUE)
      temp = temp[temp$dist!=0,]
      v_list[[paste(name, "-20_",uper_bound,"_OnlyWeakInr", sep="")]] = temp$Delta_Signal[temp$dist >-20 & temp$dist<uper_bound ]
      #v_name_list <- c(v_name_list, paste(name, "-20_0","_OnlyWeakInr", sep="") )
      
      temp=getInr2N_Dist_Delta_Signal_aroundMaxTSN(df, name, d, onlyCA = TRUE)
      temp = temp[temp$dist!=0,]
      v_list[[paste(name, "-20_",uper_bound,"_OnlyCA", sep="")]] = temp$Delta_Signal[temp$dist >-20 & temp$dist<uper_bound ]
      #v_name_list <- c(v_name_list, paste(name, "-20_0","_OnlyCA", sep="") )
      
      str(v_list)
    }
  }
}

par(mar=c(16.1, 4.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
#par(mgp=c(3,1,0))
vioplot(v_list, las=2, 
        ylim=c(-6,6),
        ylab = "log2(High Allele + 1 / Low Allele + 1)",
        col=c("black", "blue", "orange"))
abline(h=0, col="gray")



v_list_1=v_list
v_list_2=v_list
v_list_3=v_list
v_list_4=v_list

vioplot(raw_out$Delta_Signal[raw_out$dist>=-20 & raw_out$dist<0], raw_out$Delta_Signal[raw_out$dist> 0 & raw_out$dist<=20],
        names = c("BN -20-0", "BN 0-20"),      
        ylab = "log2(High Allele + 1 / Low Allele + 1)",
        names = tussue_list ,
        ylab= "log10(AT length)",
        ylim = c(0,7),
        main=Tversion )
