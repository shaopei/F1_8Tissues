setwd("~/Box Sync/KD_IGV/2020July/")
file_dir="~/Box Sync/KD_IGV/"
SNP.bw <- paste(file_dir, "P.CAST_M.B6_indelsNsnps_CAST.bam.snp.unfiltered_plus.bw", sep="")

# functions:
source("~/Box\ Sync/Danko_lab_work/F1_8Tissues/PolyA_Allele-specific/heatmap/heatmaps.R")
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

read_read_mat_SNPs <-function (SNP.bw , bed6, step=2, navg = 20, times=1, use.log=FALSE)
{
  bw.plus  <- load.bigWig(SNP.bw )
  
  hCountMatrix <- bed.step.bpQuery.bigWig(bw.plus, bed6[,c(1:3)] , step=step, abs.value=TRUE, op = "sum")
  hCountMatrix <- lapply(1:NROW(hCountMatrix), function(i){ if(bed6[i,6]=="-") return(rev(hCountMatrix[[i]])) else return(hCountMatrix[[i]])} );
  if (!use.log){
    hmat <- times * matrix(unlist(hCountMatrix), nrow= NROW(bed6), byrow=TRUE) ;
  } else {
    hmat <- log(times * matrix(unlist(hCountMatrix), nrow= NROW(bed6), byrow=TRUE) + 1) ;
  }
  #avgMat <- t(sapply(1:floor(NROW(hmat)/navg), function(x) {colMeans(hmat[((x-1)*navg+1):min(NROW(hmat),(x*navg)),])}))
  
  unload.bigWig(bw.plus);
  
  return(hmat);    
}


# pause center analysis, using dREG sites with KS test, map3 position, fdr<=0.1
combine_pause0.1 <- NULL
for (t in c("HT","SK","KD")){
  show.window=100
  pause_window_0.1 <- read.table(paste(file_dir,t,"_dREG_5mat5pat_uniq_pValue_fdr0.1.bed", sep =""), header = F)
  cat (t, dim(pause_window_0.1), "\n")
  pause_window_0.1 <- pause_window_0.1[pause_window_0.1$V1 != 'chrX',]
  cat (dim(pause_window_0.1), "\n")
  end=".rpm.bw"; times=10
  #end=".bw"; times=1
  # allelic reads (map3) #mapped position of the 3 prime end of reads
  file.bw.plus.pat <- paste(file_dir,t,".pat.map2ref.1bp_plus",end, sep="")
  file.bw.minus.pat <- paste(file_dir,t,".pat.map2ref.1bp_minus",end, sep="")
  file.bw.plus.mat <- paste(file_dir,t,".mat.map2ref.1bp_plus",end, sep="")
  file.bw.minus.mat <- paste(file_dir,t,".mat.map2ref.1bp_minus",end, sep="")
  # allelic reads (map5) HT.mat.map5.map2ref.1bp_minus.bw
  map5.file.bw.plus.pat <- paste(file_dir, "map5/",t,".pat.map5.map2ref.1bp_plus.bw", sep="")
  map5.file.bw.minus.pat <- paste(file_dir, "map5/",t,".pat.map5.map2ref.1bp_minus.bw", sep="")
  map5.file.bw.plus.mat <- paste(file_dir, "map5/",t,".mat.map5.map2ref.1bp_plus.bw", sep="")
  map5.file.bw.minus.mat <- paste(file_dir, "map5/",t,".mat.map5.map2ref.1bp_minus.bw", sep="")
  # all reads
  file.plus.bw <- paste(file_dir,t,"_PB6_F5N6_dedup_QC_end_plus",end, sep="")
  file.minus.bw <- paste(file_dir,t,"_PB6_F5N6_dedup_QC_end_minus",end, sep="")
  map5.file.plus.bw <- paste(file_dir, "map5/",t,"_PB6_F5N6_dedup_QC_end_map5_plus.bw", sep="")
  map5.file.minus.bw <- paste(file_dir, "map5/",t,"_PB6_F5N6_dedup_QC_end_map5_minus.bw", sep="")
  SNP.bw <- paste(file_dir, "P.CAST_M.B6_indelsNsnps_CAST.bam.snp.unfiltered_plus.bw", sep="")
  AT=pause_window_0.1
  # parameter setting 

  dist=200; step=1;file.pdf="heatmap.pdf"; map5=TRUE; metaplot.pdf="metaplot.pdf";
  heatmap=TRUE; metaplot=TRUE; metaplot.ylim=c(0,20);
  center.at.TSN=TRUE;
  navg = 1; use.log=FALSE; times=1
  up_dist =100;
  breaks=seq(0, 20, 0.001); times=times; map5=TRUE; metaplot.ylim = c(0,160)

  
  ### from heatmap.Pause function ###
  # AT here is the window of pause, dREG sites
  AT <- AT[,1:6] 
  
  # identify the location of max reads in mat and pat pause
  # determine the location of long pause and short pause
  hmat.pat.peak <- NULL  #map3
  hmat.mat.peak <- NULL  # map3
  map5.pat.peak <- NULL
  map5.mat.peak <- NULL
  #map5.proseq.peak <- NULL
  for (i in 1:NROW(AT)){
    # identify the base with of max map3 reads 
    hmat.pat <- read_read_mat_S (file.bw.plus.pat, file.bw.minus.pat, AT[i,], step, times=times, use.log=use.log)
    hmat.mat <- read_read_mat_S (file.bw.plus.mat, file.bw.minus.mat, AT[i,], step, times=times, use.log=use.log) 
    hmat.pat.peak[i] <- which(hmat.pat == max(hmat.pat))
    hmat.mat.peak[i] <- which(hmat.mat == max(hmat.mat))
    # identify the base with of max map5 reads 
    map5.pat <- read_read_mat_S (map5.file.bw.plus.pat, map5.file.bw.minus.pat, AT[i,], step, times=times, use.log=use.log)
    map5.mat <- read_read_mat_S (map5.file.bw.plus.mat, map5.file.bw.minus.mat, AT[i,], step, times=times, use.log=use.log) 
    map5.pat.peak[i] <- which(map5.pat == max(map5.pat))  # if two bases with same read count, output the 5 prime one (up stream)
    map5.mat.peak[i] <- which(map5.mat == max(map5.mat)) 
    
    map5.t=c(map5.pat.peak[i],   map5.mat.peak[i])   # the base with of max map5 reads 
    
    # determine which allele is early pause, which is late pause
    a=sort.int(c(hmat.pat.peak[i],  hmat.mat.peak[i]), index.return = TRUE )
    AT$early.pause[i] = a$x[1] 
    AT$late.pause[i] = a$x[2]
    AT$TSN.early.pause[i] = map5.t[a$ix[1]]  # from the same allele as early pause, which base has the max map5 reads
    AT$TSN.late.pause[i] = map5.t[a$ix[2]]   # from the same allele as late pause, which base has the max map5 reads,
  }
  AT$tissue=t
  ### from heatmap.Pause function end ###
  combine_pause0.1=rbind.data.frame(combine_pause0.1, AT)
}

AT=combine_pause0.1

# distribution of distance between early and late TSN
AT$early.TSN.pause.dist = AT$early.pause - AT$TSN.early.pause
AT$late.TSN.pause.dist = AT$late.pause - AT$TSN.late.pause

dim(AT)

TSN.pause.dist.lowerbound=10
TSN.pause.dist.upperbpund=50

AT$valid.pairs = AT$early.TSN.pause.dist > TSN.pause.dist.lowerbound & AT$late.TSN.pause.dist > TSN.pause.dist.lowerbound &
  AT$early.TSN.pause.dist <TSN.pause.dist.upperbpund & AT$late.TSN.pause.dist <TSN.pause.dist.upperbpund
sum(AT$valid.pairs)

subAT=AT[!AT$valid.pairs,]
AT = AT[AT$valid.pairs,]

AT$TSN.dist = AT$TSN.late.pause - AT$TSN.early.pause
#View(AT)
hist(AT$TSN.dist, main="distance between early and late TSN")
hist(AT$TSN.dist[AT$TSN.dist !=0 ], main="distance between early and late TSN, exclude dist==0")
sum(AT$TSN.dist == 0)
sum(AT$TSN.dist !=0 )
sum(AT$TSN.dist <0 )
sum(AT$TSN.dist >0 )
sum(AT$TSN.dist >10 )

# distribution of distance between early and late pause
AT$pause.dist= AT$late.pause - AT$early.pause
hist(AT$pause.dist, main="distance between early and late pause")
hist(AT$pause.dist[AT$pause.dist!= 0], main="distance between early and late pause, exclude dist==0",
     add=T, col="red")
sum(AT$pause.dist==0)
sum(AT$pause.dist>0)
sum(AT$pause.dist> 10)

#hist(AT$pause.dist[AT$TSN.dist== 0], add=T, col="blue", density = 45, breaks=seq(0,400,50))


hist(AT$pause.dist - AT$TSN.dist
     , breaks=seq(-35.5,35,1)
)
#panel52_TSN_Pause_hist_plot.pdf
# Figure 4B
pdf("TSN_Pause_hist.pdf", width=7, height = 7)
par(mfrow=c(2,1))
par(mar=c(6.1, 7.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
par(mgp=c(3,1,0))
par(cex.lab=2.2, cex.axis=2.2)
hist(abs(AT$TSN.dist[AT$pause.dist==0])# &AT$TSN.dist !=0 ]
     , breaks=seq(-0.5,35,1), col="gray", las=1
     #, ylim=c(0,20)
)
hist(AT$pause.dist[AT$TSN.dist==0]
     , breaks=seq(-0.5,35,1), col="gray", las=1
)
dev.off()
sum(AT$pause.dist==0)
sum(AT$pause.dist==0 &AT$TSN.dist !=0)

sum(AT$TSN.dist==0)
sum(AT$TSN.dist!=0)


subAT=AT[abs(AT$pause.dist - AT$TSN.dist) > 10,]
subAT=subAT[subAT$TSN.dist==0,]
subAT$deltaTsnPauseDist = subAT$pause.dist - subAT$TSN.dist

library("ggplot2")
# Figure 4A
#panel51_TSN_Pause_scatterplot.pdf
pdf("TSN_Pause_scatterplot.pdf", width=7, height = 7)
myColor <- rev(RColorBrewer::brewer.pal(11, "Spectral"))
myColor_scale_fill <- scale_fill_gradientn(colours = myColor, trans='log10')

p <- ggplot(AT, aes(x=TSN.dist, y=pause.dist)) +
  geom_bin2d(bins = 100) +  myColor_scale_fill + theme_bw() +  
  #xlim(-20, 50) + ylim(0,50) + #labs(x="log10 (TSS Read Counts)")  + labs(y="TSN counts") +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=20), #,face="bold"))
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 14)
  )
p
dev.off()

### maxTSN centered pause analysis ###
# use the shared allelic maxTSN
# indel analysis
#Tissue="HT"
df<-NULL
for (Tissue in c("HT","SK", "KD")){
  df1=read.table(paste(Tissue, "_BothAlleleMaxTSNs_ratio0.5-2_RLmatpat_map3refTomaxTSNmatpat_ClosetIndel.bed", sep = ""))
  df1$Tissue=Tissue
  df = rbind.data.frame(df,df1)
}

#colnames(df)[7:10]=c("mat_RL", "pat_RL" , "mat_map3RefDist", "pat_map3RefDist" )
#View(df)
for (i in 1:dim(df)[1]){
  # pick maxPause
  # if there is a tie, the shorter read length is reported
  df$mat_maxPause_RL[i] = as.numeric(names(sort(table(unlist(strsplit(as.character(df$V7[i]), ","))),decreasing=TRUE)[1]))
  df$pat_maxPause_RL[i] = as.numeric(names(sort(table(unlist(strsplit(as.character(df$V8[i]), ","))),decreasing=TRUE)[1]))
  df$mat_maxPause_map3[i] = as.numeric(names(sort(table(unlist(strsplit(as.character(df$V9[i]), ","))),decreasing=TRUE)[1]))
  df$pat_maxPause_map3[i] = as.numeric(names(sort(table(unlist(strsplit(as.character(df$V10[i]), ","))),decreasing=TRUE)[1]))
  df$mat_Idel_length[i]= length(unlist(strsplit(unlist(strsplit(as.character(df$V14[i]), ","))[1], "")))
  df$pat_Idel_length[i]= length(unlist(strsplit(unlist(strsplit(as.character(df$V14[i]), ","))[2], "")))
  df$maxPauseSite_map3[i] = max(df$mat_maxPause_map3[i],  df$pat_maxPause_map3[i])
  df$maxPauseSite_RL[i] = max(df$mat_maxPause_RL[i],  df$pat_maxPause_RL[i])
  
  df$earlyPause[i]= min(df$mat_maxPause_map3[i],df$pat_maxPause_map3[i])
  df$latePause[i]= max(df$mat_maxPause_map3[i],df$pat_maxPause_map3[i])
  
  }



df$Indel_Len = df$mat_Idel_length - df$pat_Idel_length

# KS test of pause shape, using read length RL (actual length on native genome) and map3 mapped position on reference genome
for (i in 1:dim(df)[1]){
  m=as.numeric(unlist(strsplit(as.character(df$V7[i]), ","))) #mat read length, between maxTSN and pause
  p=as.numeric(unlist(strsplit(as.character(df$V8[i]), ","))) # pat read length
  df$RL.p.value[i] = ks.test(m,p) $ p.value
}
for (i in 1:dim(df)[1]){
  m=as.numeric(unlist(strsplit(as.character(df$V9[i]), ","))) # mat reads 3 prime mapped position on mm10
  p=as.numeric(unlist(strsplit(as.character(df$V10[i]), ","))) # pat reads 3 prime mapped position on mm10
  df$map3.p.value[i] = ks.test(m,p) $ p.value
}

df$RL.p.value.fdr = p.adjust(df$RL.p.value, method = "fdr")
df$map3.p.value.fdr = p.adjust(df$map3.p.value, method = "fdr")

# indel length VS map3 position difference
# V15 is the distance of the indel to maxTSN
lim=40

## average
for (i in 1:dim(df)[1]){
  df$mat_AvePause_RL[i] = mean(as.numeric(unlist(strsplit(as.character(df$V7[i]), ","))))
  df$pat_AvePause_RL[i] = mean(as.numeric(unlist(strsplit(as.character(df$V8[i]), ","))))
  df$mat_AvePause_map3[i] = mean(as.numeric(unlist(strsplit(as.character(df$V9[i]), ","))))
  df$pat_AvePause_map3[i] = mean(as.numeric(unlist(strsplit(as.character(df$V10[i]), ","))))
}

# df$V15 is the distance between maxTSN to the nearest indel
# df$V15 <= df$maxPauseSite_map3  to identify the indels located between maxPauseSite_map3 and maxTSN
# (df$mat_maxPause_map3 != df$pat_maxPause_map3) make sure the early pause position =! late pause position
# df$map3.p.value.fdr <=0.1 KS test is significant (fdr<=0.1)

# use only fdr<=0.1 sites where early pause != late pause
df$target = (df$V15 <= df$maxPauseSite_map3) & (df$map3.p.value.fdr <=0.1 ) & (df$mat_maxPause_map3 != df$pat_maxPause_map3)
dim(df)
sum(df$target)
#figure 4C
# Panel53_pause_indelength.pdf
pdf("pause_indelength.pdf", width=7, height = 7, useDingbats=FALSE)
par(mar=c(6.1, 7.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
par(mgp=c(3,1,0))
par(cex.lab=2.2, cex.axis=2.2)
pch_u=19
plot( (df$mat_AvePause_map3- df$pat_AvePause_map3)[df$target& df$Tissue=="HT"],
      (df$mat_Idel_length - df$pat_Idel_length)[df$target & df$Tissue=="HT"],
      xlim=c(-1*lim,lim),
      ylim=c(-1*lim,lim),
      xlab="B6 - CAST average pause position on mm10",
      ylab="B6 - CAST indel length",
      pch=pch_u, frame=F,
      cex=2.2, las=1,
      col="red")
# points( (df$mat_AvePause_map3- df$pat_AvePause_map3)[df$target& df$Tissue=="HT"], 
#         (df$mat_Idel_length - df$pat_Idel_length)[df$target & df$Tissue=="HT"],
#         pch=pch_u, 
#         col="red")
points( (df$mat_AvePause_map3- df$pat_AvePause_map3)[df$target& df$Tissue=="SK"], 
        (df$mat_Idel_length - df$pat_Idel_length)[df$target & df$Tissue=="SK"],
        pch=pch_u, 
        cex=2.2,
        col="dark orange")
points( (df$mat_AvePause_map3- df$pat_AvePause_map3)[df$target& df$Tissue=="KD"], 
        (df$mat_Idel_length - df$pat_Idel_length)[df$target & df$Tissue=="KD"],
        pch=pch_u,
        cex=2.2, 
        col="blue")

abline(a=5, b=-1, col="gray")
abline(a=0,b=-1, col="dark gray")
abline(a=-5, b=-1, col="gray")

legend("topright", 
       legend = c("HT", "SK", "KD"),
       pch=pch_u,
       #       cex=2, 
       lty=c(0,0),
       #bty="n",
       #lwd=1.5, 
       #density=c(10000,25),
       #angle=c(180,45),
       #angle=45,
       cex=2.2,
       col=c("red", "dark orange", "blue")
       , bty = "n"
       
)
dev.off()

# return the row with distinct chromStart and early Pause sites
# UniqRowCount <- function(subdf){
#   dim(unique(data.frame(subdf$V2, subdf$earlyPause)))[1]
# }
# return the number of rows with uniq early pause position
UniqRowCount <- function(subdf){
  bed6 <- subdf[1:6]
  if (NROW(bed6)>0){
  for (i in 1:NROW(bed6)){
    if(bed6[i,6]=="-") {
      bed6[i,3] <- subdf[i,3] - subdf$earlyPause[i]
      bed6[i,2] <- subdf[i,2] - subdf$earlyPause[i]
    } else {
      bed6[i,2] <- subdf[i,2] + subdf$earlyPause[i]
      bed6[i,3] <- subdf[i,3] + subdf$earlyPause[i]
    }
  }
  
  bed6[,4]=111
  bed6[,5]=111
  dim(unique(bed6))[1]
  } else {0}
  #dim(unique(data.frame(subdf$V1, subdf$V2+subdf$earlyPause)))[1]
}


# Fisher's exact test, test if there is more sites contains indel in the group of   ks test fdr<0.1 AND early pause != late pause
sum((df$map3.p.value.fdr <=0.1 ) & (df$mat_maxPause_map3 != df$pat_maxPause_map3)) #308
a=UniqRowCount(df[(df$map3.p.value.fdr <=0.1 ) & (df$mat_maxPause_map3 != df$pat_maxPause_map3),])  #use maxTSN and early pause postion to count uniq row
b=UniqRowCount(df[df$V15 <= df$maxPauseSite_map3 & (df$map3.p.value.fdr <=0.1 ) & (df$mat_maxPause_map3 != df$pat_maxPause_map3),])
c=UniqRowCount(df[(df$map3.p.value.fdr >0.9 ) & (df$mat_maxPause_map3 == df$pat_maxPause_map3),])
d=UniqRowCount(df[df$V15 <= df$maxPauseSite_map3 &(df$map3.p.value.fdr >0.9 ) & (df$mat_maxPause_map3 == df$pat_maxPause_map3),])
a;b;c;d
#269;56;1396;113
#                                                  # of sites with indel  VS # of sites 
# ks test fdr<0.1 AND early pause != late pause       56                  269
# ks test fdr>0.9 AND early pause == late pause      113                 1396   
fisher.test(data.frame(c(b, a), c(d, c)))
#p-value = 0.0000003912

### maxTSN centered pause analysis ###
# use the shared allelic maxTSN
# SNP analysis
df1=read.table("HT_BothAlleleMaxTSNs_ratio0.5-2_map3_mat.pat.IDEreads_map2ref_DistanceTomaxTSN.bed")
df2=read.table("SK_BothAlleleMaxTSNs_ratio0.5-2_map3_mat.pat.IDEreads_map2ref_DistanceTomaxTSN.bed")
df3=read.table("KD_BothAlleleMaxTSNs_ratio0.5-2_map3_mat.pat.IDEreads_map2ref_DistanceTomaxTSN.bed")

df1$Tissue="HT"
df2$Tissue="SK"
df3$Tissue="KD"
df=rbind.data.frame(df1,df2,df3)  

colnames(df)[7:9]=c("mat_map3", "pat_map3" , "ide_map3" )
#View(df)
for (i in 1:dim(df)[1]){
  # pick maxPause
  # if there is a tie, the shorter read length/mapping position is reported
  reads=as.numeric(unlist(c(strsplit(as.character(df$mat_map3[i]), ","), strsplit(as.character(df$pat_map3[i]), ","), strsplit(as.character(df$ide_map3[i]), ","))))
  df$maxPause_map3AllReads[i] = as.numeric(names(sort(table(reads),decreasing=TRUE)[1]))
  if (!is.na(df$mat_map3[i])){
    df$mat_maxPause_map3[i] = as.numeric(names(sort(table(unlist(strsplit(as.character(df$mat_map3[i]), ","))),decreasing=TRUE)[1]))
  }
  if (!is.na(df$pat_map3[i])){
    df$pat_maxPause_map3[i] = as.numeric(names(sort(table(unlist(strsplit(as.character(df$pat_map3[i]), ","))),decreasing=TRUE)[1]))
    df$earlyPause[i]= min(df$mat_maxPause_map3[i],df$pat_maxPause_map3[i])
    df$latePause[i]= max(df$mat_maxPause_map3[i],df$pat_maxPause_map3[i])
  }
  
}

# with or without indel 
Tissue_list= c("HT", "SK", "KD")
indel<- NULL
for (i in 1:3){
  indel_i=read.table(paste(Tissue_list[i], "_BothAlleleMaxTSNs_ratio0.5-2_RLmatpat_map3refTomaxTSNmatpat_ClosetIndel.bed", sep = ""))
  indel_i$Tissue=Tissue_list[i]
  indel=rbind.data.frame(indel, indel_i)  
}

indel$V5=111
colnames(indel)[14:15]=c("indel", "dist_Indel_maxTSN")
df_indel <- merge(df, indel, sort=F)
df = df_indel


dim(df)
sum(df$Tissue=="HT")
sum(df$Tissue=="SK")
sum(df$Tissue=="KD")

# get the position of SNPs around the early pause
show.window=20; step=1
bed6 <- df[,1:6]
for (i in 1:NROW(bed6)){
  if(bed6[i,6]=="-") {
    bed6[i,3] <- df[i,3] - df$earlyPause[i] + show.window
    bed6[i,2] <- df[i,2] - df$earlyPause[i] - show.window
  } else {
    bed6[i,2] <- df[i,2] + df$earlyPause[i] - show.window
    bed6[i,3] <- df[i,3] + df$earlyPause[i] + show.window
  }
}


#write.table(bed6, file="test.bed", quote = F, sep="\t", row.names = F, col.names = F)
SNPs.show.window <- read_read_mat_SNPs (SNP.bw, bed6[,c(1:6)], step=step, times=1, use.log=FALSE)
dim(SNPs.show.window)
for (i in 1:NROW(bed6)){
  SNP_distance_earlyPause = which(SNPs.show.window[i,]==1)-show.window-1
  df$SNP_distance_earlyPause[i] = paste(SNP_distance_earlyPause,collapse = ",")
  df$SNP1_distance_earlyPause[i] = SNP_distance_earlyPause[order(abs(SNP_distance_earlyPause), decreasing = F)[1]]
  df$SNP2_distance_earlyPause[i] = SNP_distance_earlyPause[order(abs(SNP_distance_earlyPause), decreasing = F)[2]]
  df$SNP3_distance_earlyPause[i] = SNP_distance_earlyPause[order(abs(SNP_distance_earlyPause), decreasing = F)[3]]
}


# KS test result
df$map3.p.value=1
for (i in 1:dim(df)[1]){
  m=as.numeric(unlist(strsplit(as.character(df$mat_map3[i]), ",")))
  p=as.numeric(unlist(strsplit(as.character(df$pat_map3[i]), ",")))
  if (!is.na(m) & !is.na(p) ){
    df$map3.p.value[i] = ks.test(m,p) $ p.value
  }
}

df$map3.p.value.fdr = p.adjust(df$map3.p.value, method = "fdr")

df$earlyPauseAllele = "CAST"
df$earlyPauseAllele[df$earlyPause==df$mat_maxPause_map3] = "CAST"

write.table(df[,c("V1","V2", "V3", "V4", "V5", "V6", "Tissue" ,"earlyPause", "latePause", "map3.p.value.fdr", "earlyPauseAllele")], file="T3_BothAlleleMaxTSNs_ratio0.5-2_map3_mat.pat.IDEreads_map2ref_DistanceTomaxTSN.bed", quote = F, sep="\t",
            row.names = F, col.names = T)

getwd()

metaplot.SNPsLocation.aroundMaxPause <-function(df, name="", use.sum=FALSE, col="red", show.window = 49, step=1 ,add=FALSE, pch_u=19){
  bed6 <- df[,1:6]
  # maxPause location +- window
  
  for (i in 1:NROW(bed6)){
    if(bed6[i,6]=="-") {
      bed6[i,3] <- df[i,3] - df$maxPause_map3AllReads[i] + show.window
      bed6[i,2] <- df[i,2] - df$maxPause_map3AllReads[i] - show.window
    } else {
      bed6[i,2] <- df[i,2] + df$maxPause_map3AllReads[i] - show.window
      bed6[i,3] <- df[i,3] + df$maxPause_map3AllReads[i] + show.window
    }
  }
  
  bed6[,4]=0
  bed6[,5]=111
  bed6 = unique(bed6) # remove duplicates with max pause position  
  
  #write.table(bed6, file="test.bed", quote = F, sep="\t", row.names = F, col.names = F)
  SNPs.show.window <- read_read_mat_SNPs (SNP.bw, bed6[,c(1:6)], step=step, times=1, use.log=FALSE)
  
  
  if (use.sum){
    a = colSums(SNPs.show.window )
  }else{
    a = colMeans(SNPs.show.window )
  }
  
  
  #pdf(metaplot.pdf, width=6, height = 6)
  #par(mar=c(6.1, 7.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
  #par(mgp=c(3,1,0))
  #par(cex.lab=2.2, cex.axis=2.2)
  
  if (step >1){
    y <- a
    x <- seq(-1*(show.window)+step/2, show.window, step)[1:length(y)]
  }else{
    y <- a
    x <- seq(-1*(show.window), show.window, step)[1:length(y)]
  }
  if(add){
    lines(x, y, col=col, type="o", pch=pch_u)
  }else{
    plot(x, y, col=col, 
         main = name,
         ylab= "SNPs counts mean",
         type = "o", 
         pch=pch_u,
         xlab="",
         las=1
         #ylim = c(0,0.5)
         
    )
  }
  
  #dev.off()
  return (list(x, y))
}


metaplot.SNPsLocation.aroundEarlyPause <-function(df, name="", use.sum=FALSE, col="red", show.window = 49, step=1 ,add=FALSE, plot=TRUE, pch_u=19){
  bed6 <- df[,1:6]
  bed6[,4] <-df$earlyPause
  #dim(bed6)
  bed6_copy = bed6
  # maxPause location +- window
  
  for (i in 1:NROW(bed6)){
    if(bed6[i,6]=="-") {
      bed6[i,3] <- bed6_copy[i,3] - bed6_copy[i,4] + show.window
      bed6[i,2] <- bed6_copy[i,2] - bed6_copy[i,4] - show.window
    } else {
      bed6[i,2] <- bed6_copy[i,2] + bed6_copy[i,4] - show.window
      bed6[i,3] <- bed6_copy[i,3] + bed6_copy[i,4] + show.window
    }
  }
  
  bed6[,4]=0
  bed6[,5]=111
  bed6 = unique(bed6) # remove duplicates with max pause position  
  
  #write.table(bed6, file="test.bed", quote = F, sep="\t", row.names = F, col.names = F)
  SNPs.show.window <- read_read_mat_SNPs (SNP.bw, bed6[,c(1:6)], step=step, times=1, use.log=FALSE)
  dim(SNPs.show.window)
 # for (i in 1:NROW(bed6)){
#    SNP_distance_earlyPause = which(SNPs.show.window[i,]==1)-show.window-1
    
 #   df$SNP1_distance_earlyPause[i] = SNP_distance_earlyPause[order(abs(SNP_distance_earlyPause), decreasing = F)[1]]
#    df$SNP2_distance_earlyPause[i] = SNP_distance_earlyPause[order(abs(SNP_distance_earlyPause), decreasing = F)[2]]
#    df$SNP3_distance_earlyPause[i] = SNP_distance_earlyPause[order(abs(SNP_distance_earlyPause), decreasing = F)[3]]
#  }
  
  
  if (use.sum){
    a = colSums(SNPs.show.window )
  }else{
    a = colMeans(SNPs.show.window )
  }
  
  
  #pdf(metaplot.pdf, width=6, height = 6)
  #par(mar=c(6.1, 7.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
  #par(mgp=c(3,1,0))
  #par(cex.lab=2.2, cex.axis=2.2)
  
  if (step >1){
    y <- a
    x <- seq(-1*(show.window)+step/2, show.window, step)[1:length(y)]
  }else{
    y <- a
    x <- seq(-1*(show.window), show.window, step)[1:length(y)]
  }
  if(plot){
  if(add){
    lines(x, y, col=col, type="o", pch=pch_u)
  }else{
    plot(x, y, col=col, 
         main = name,
         ylab= "SNPs counts mean",
         type = "o", 
         xlab="distance to early pause",
         pch=pch_u,
         las=1
         #ylim = c(0,0.5)
         
    )
  }
  }
  
  #dev.off()
  return (list(x, y, dim(bed6)[1]))
}


metaplot.SNPsLocation.aroundLatePause <-function(df, name="", use.sum=FALSE, col="red", show.window = 49, step=1 ,add=FALSE,pch_u=19){
  bed6 <- df[,1:6]
  # maxPause location +- window
  
  for (i in 1:NROW(bed6)){
    if(bed6[i,6]=="-") {
      bed6[i,3] <- df[i,3] - df$latePause[i] + show.window
      bed6[i,2] <- df[i,2] - df$latePause[i] - show.window
    } else {
      bed6[i,2] <- df[i,2] + df$latePause[i] - show.window
      bed6[i,3] <- df[i,3] + df$latePause[i] + show.window
    }
  }
  
  bed6[,4]=0
  bed6[,5]=111
  bed6 = unique(bed6) # remove duplicates with max pause position  
  
  #write.table(bed6, file="test.bed", quote = F, sep="\t", row.names = F, col.names = F)
  SNPs.show.window <- read_read_mat_SNPs (SNP.bw, bed6[,c(1:6)], step=step, times=1, use.log=FALSE)
  
  
  if (use.sum){
    a = colSums(SNPs.show.window )
  }else{
    a = colMeans(SNPs.show.window )
  }
  
  
  #pdf(metaplot.pdf, width=6, height = 6)
  #par(mar=c(6.1, 7.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
  #par(mgp=c(3,1,0))
  #par(cex.lab=2.2, cex.axis=2.2)
  
  if (step >1){
    y <- a
    x <- seq(-1*(show.window)+step/2, show.window, step)[1:length(y)]
  }else{
    y <- a
    x <- seq(-1*(show.window), show.window, step)[1:length(y)]
  }
  if(add){
    lines(x, y, col=col, type="o", pch=pch_u)
  }else{
    plot(x, y, col=col, 
         main = name,
         ylab= "SNPs counts mean",
         type = "o", 
         xlab="",
         pch=pch_u,
         las=1
         #ylim = c(0,0.5)
         
    )
  }
  
  #dev.off()
  return (list(x, y))
}



### early pause with at least bp.apart bp difference
# exclude region with indel
tempFunc <-function(bp.apart=1, upto=100, show.window=30){
  step=1
  
  ## SNPs around early pause
  t1=metaplot.SNPsLocation.aroundEarlyPause(df[df$map3.p.value.fdr<=0.1 
                                            & abs(df$mat_maxPause_map3 - df$pat_maxPause_map3) >= bp.apart 
                                            & abs(df$mat_maxPause_map3 - df$pat_maxPause_map3) <= upto,], 
                                            #& df$dist_Indel_maxTSN > df$latePause,], 
                                            show.window=show.window,
                                         name=paste("HSK, step=", step,", early and late at least ", bp.apart," bp apart, upto " ,upto, sep=""), col="red", step=step)
  abline(v=0, col="gray")

  t2=metaplot.SNPsLocation.aroundEarlyPause(df[df$map3.p.value.fdr>0.9  
                                               & (df$mat_maxPause_map3 == df$pat_maxPause_map3) ,], 
                                            #& df$dist_Indel_maxTSN > df$latePause,], 
                                            show.window=show.window,
                                         col="blue", step=step, add=TRUE)

  legend("topright", legend=c(paste("Allelic pause difference, n= ", t1[[3]], sep=""),
                              paste("No allelic pause difference, n= ", t2[[3]], sep="")),
         col=c("red", "blue"), 
         bty = "n", lty=1, pch=19)

}

##
## SNPs around early pause, use.sum = TRUE, for test
step=1; bp.apart=1; upto=100
e1=metaplot.SNPsLocation.aroundEarlyPause(df[df$map3.p.value.fdr<=0.1 
                                             & abs(df$mat_maxPause_map3 - df$pat_maxPause_map3) >= bp.apart 
                                             & abs(df$mat_maxPause_map3 - df$pat_maxPause_map3) <= upto,], 
                                             #& df$dist_Indel_maxTSN > df$latePause,], 
                                          name=paste("HSK, step=", step,", early and late at least ", bp.apart," bp apart, upto " ,upto, sep=""), col="red", step=step,
                                          use.sum = TRUE)
abline(v=0, col="gray")
# abline(v=-5, col="gray")
# abline(v=-10, col="gray")

e9=metaplot.SNPsLocation.aroundEarlyPause(df[df$map3.p.value.fdr>0.9
                                             & (df$mat_maxPause_map3 == df$pat_maxPause_map3) ,], 
                                             #& df$dist_Indel_maxTSN > df$latePause,], 
                                            col="blue", step=step, use.sum = TRUE)
# legend("topright", legend=c(paste("fdr <= 0.1, n= ", sum(df$map3.p.value.fdr<=0.1 & abs(df$mat_maxPause_map3 - df$pat_maxPause_map3) >0), sep=""), 
#                             paste("fdr >  0.9, n= ", sum(df$map3.p.value.fdr>0.9), sep="")),
#        col=c("red", "blue"), bty = "n", lty=1, pch=19)
# 

# add a test
#pdf(paste(organ,"_FisherExactTest_CountOfTSSwithAT2GCSNPs_vs_AllTSS.pdf", sep=""))
# panel54_SNP_distribution_around_short_pause
#Figure4D
pdf("SNP_distribution_around_short_pause.pdf", width=9.25, height = 5.36, useDingbats=FALSE)
par(mar=c(6.1, 7.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
par(mgp=c(3,1,0))
par(cex.lab=2.2, cex.axis=2.2)
tempFunc(1,100, show.window=30) # sites with indel were removed
dev.off()
par(mfrow=c(3,1))
tempFunc(1,100) # use mean
fdr_cutoff=0.1

Test_one_by_one=FALSE
if(Test_one_by_one){
  p.value <- NULL
  odds.ratio <- NULL
  testors <- NULL
  
  for (n in (1:length(e1[[2]]))){
    
    testor <-  matrix(c(e1[[2]][n], e1[[3]],
                        e9[[2]][n], e9[[3]]),
                      nrow = 2,
                      dimnames = list(TSS = c("with SNPs", "with or without SNPs"),
                                      KS.Test = c("fdr <= 0.1, Early Pause", "fdr >0.9, Early Pause "))); testor
    
    f = fisher.test(testor, alternative = "greater" ); f
    p.value= c(p.value, f$p.value)
    odds.ratio = c(odds.ratio, f$estimate)
  }
  adjust.p = p.adjust(p.value, method="fdr")
  plot(e1[[1]], -log10(adjust.p), type="o" , xlab="distance to early pause", main=" pause KS test fdr <= 0.1 vs fdr > 0.9", las=1, col="red")
  abline(h=-1*log10(fdr_cutoff),col="gray")
  plot(e1[[1]],odds.ratio, type="o" , xlab="distance to early pause", las=1)
  abline(h=1,col="gray")
  #plot(e1[[1]], adjust.p, type="o" , #ylim = c(0,0.2),
  #     xlab="distance to early pause", main=" pause KS test fdr <= 0.1 vs fdr > 0.9", las=1, col="red")
  #abline(h=0.25,col="gray")
  # text(e1[[1]][which(p.value <= 0.1)],p.value[which(p.value <= 0.1)], label=paste(round(p.value[which(p.value <= 0.1)], digits = 2), sep=" "))
  # dev.off()
}

group_fisher_test <- function(group1,e1,e9){
  e1_i=0 ; e9_i=0
  for (i in which(e1[[1]]%in%group1)){
    e1_i = e1_i + e1[[2]][i]
    e9_i = e9_i + e9[[2]][i]
  }
  return(fisher.test(data.frame(c(e1_i, e1[[3]]),c(e9_i, e9[[3]])), alternative = "greater" ))
}

if(!Test_one_by_one){
  group1_list=list(c(-4:-8), c(-2,-3), c(0), c(1))#, c(2:7), 8)
  group1_label=c("-4-8", "-2,-3", 0, 1)#, "2:7", 8)
  #group1_label=c((-4-8)/2, (-2-3)/2,0,1, (2+10)/2)
  group1_pvalue<-NULL
  p.value <- NULL
  odds.ratio <- NULL
  for (group1 in group1_list){
    f=group_fisher_test(group1,e1,e9)
    p.value= c(p.value, f$p.value)
    odds.ratio = c(odds.ratio, f$estimate)
  }
  adjust.p = p.adjust(p.value, method="fdr")
  barplot( -log10(adjust.p), names.arg = group1_label, 
           ylab="-log10(adjust.p)" , xlab="distance to early pause", main=" pause KS test fdr <= 0.1 vs fdr > 0.9", las=1, col="red")
  abline(h=-1*log10(fdr_cutoff),col="gray")
  barplot(odds.ratio, names.arg = group1_label,
          ylab="odds.ratio" , xlab="distance to early pause", las=1)
  abline(h=1)
}

# effect of SNPs at the difference between maxPauses of early and late pause
# exclude allelic maxPause difference 0
df$AllelicMaxPauseDist = abs(df$mat_maxPause_map3 - df$pat_maxPause_map3)
# Figure 4H
sub_df=df[df$map3.p.value.fdr<=0.1 & df$AllelicMaxPauseDist >0,]
sub_df_SNP=df[df$map3.p.value.fdr<=0.1 & df$AllelicMaxPauseDist >0 & 
                df$SNP1_distance_earlyPause==0 & !is.na(df$SNP1_distance_earlyPause),]
#UniqRowCount(sub_df)
#UniqRowCount(sub_df_SNP)
dim(sub_df)
dim(sub_df_SNP)

useOnlySitesWithSNP_C2ATG=TRUE
if(useOnlySitesWithSNP_C2ATG){
  bed6_SNP <- sub_df_SNP[,c('V1','V2','V3','earlyPause','V5','V6','AllelicMaxPauseDist','V1','V2','V3')]
  bed6_SNP = unique(bed6_SNP)  # remove duplicates with shared maxTSN, early pause, and AllelicMaxPauseDist
  bed6_SNP_copy = bed6_SNP
  # maxPause location +- window
  show.window=0
  for (i in 1:NROW(bed6_SNP)){
    if(bed6_SNP[i,6]=="-") {
      bed6_SNP[i,3] <- bed6_SNP_copy[i,3] - bed6_SNP_copy[i,4] + show.window
      bed6_SNP[i,2] <- bed6_SNP_copy[i,2] - bed6_SNP_copy[i,4] - show.window
    } else {
      bed6_SNP[i,2] <- bed6_SNP_copy[i,2] + bed6_SNP_copy[i,4] - show.window
      bed6_SNP[i,3] <- bed6_SNP_copy[i,3] + bed6_SNP_copy[i,4] + show.window
    }
  }
  write.table(bed6, file="sub_df_SNP.bed", quote = F, sep="\t", row.names = F, col.names = F)
  # find out the SNP at the early pause using shell script in F1_TSN_identifyTSS_SingleBaseRunOn_forManuscript.sh
}
dim(bed6_SNP)

remove_duplicates_V2_earlyPause <- function(subdf){
  return(subdf[!duplicated(subdf[,c('V1','V2','earlyPause','V6')]),])
}

getEarlyPauseSites <- function(subdf){
  bed6 <- subdf[1:6]
  for (i in 1:NROW(bed6)){
    if(bed6[i,6]=="-") {
      bed6[i,3] <- subdf[i,3] - subdf$earlyPause[i]
      bed6[i,2] <- subdf[i,2] - subdf$earlyPause[i]
    } else {
      bed6[i,2] <- subdf[i,2] + subdf$earlyPause[i]
      bed6[i,3] <- subdf[i,3] + subdf$earlyPause[i]
    }
  }
  bed6[,4]=111
  bed6[,5]=111
  return((bed6))
  #dim(unique(data.frame(subdf$V1, subdf$V2+subdf$earlyPause)))[1]
}

set_remove_duplicates_V2_earlyPause =FALSE
if (set_remove_duplicates_V2_earlyPause){
  sub_df=remove_duplicates_V2_earlyPause(sub_df)
  sub_df_SNP=remove_duplicates_V2_earlyPause(sub_df_SNP)
}else{
  # deduplicates based on early pause sites and AllelicMaxPauseDist
  sub_df$earlyPauseSites=getEarlyPauseSites(sub_df)$V2
  sub_df=sub_df[!duplicated(sub_df[,c('V1','earlyPauseSites','V6','AllelicMaxPauseDist')]),]
  sub_df_SNP$earlyPauseSites=getEarlyPauseSites(sub_df_SNP)$V2
  sub_df_SNP=sub_df_SNP[!duplicated(sub_df_SNP[,c('V1','earlyPauseSites','V6','AllelicMaxPauseDist')]),]
}
dim(sub_df)
dim(sub_df_SNP)
# if only use those with C to ATG SNPs
if(useOnlySitesWithSNP_C2ATG){
  sub_df_SNP_2=read.table("sub_df_SNP_C2ATG.bed", header=T)
  sub_df_SNP=cbind.data.frame(bed6_SNP, sub_df_SNP_2)
  dim(sub_df_SNP)
  sub_df_SNP=sub_df_SNP[sub_df_SNP$V2==sub_df_SNP$chrStart,] # examine if the rows in bed6_SNP and sub_df_SNP_2 matches
  dim(sub_df_SNP)
  View(sub_df_SNP)
  sub_df_SNP=sub_df_SNP[sub_df_SNP$C2ATG_SNPs==1,] # if only use those with C to ATG SNPs
  dim(sub_df_SNP)
  # remove duplicates based on early pause sites and AllelicMaxPauseDist 
  sub_df_SNP = sub_df_SNP[!duplicated(sub_df_SNP[,c('V1','V6','earlyPauseV2','AllelicMaxPauseDist')]),]
  dim(sub_df_SNP)
}
#par(mfrow=c(2,1))
# panel55_SNPc_allelicMaxPauseDist
pdf("SNPc_allelicMaxPauseDist.pdf", width=7, height = 3.5, useDingbats=FALSE)
par(mar=c(6.1, 7.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
par(mgp=c(3,1,0))
par(cex.lab=2.2, cex.axis=2.2)
plot(ecdf(sub_df$AllelicMaxPauseDist), col="blue", 
     xlab="Distance between allelic maxPause",
     ylab="fraction",
     las=1, 
     main="maxPause map3TomaxTSNs KS FDR <= 0.1, e and l >1bp move"
     )
hist(sub_df$AllelicMaxPauseDist,
     freq = F, ylim=c(0,0.4), las=1,
     breaks = seq(-0.5,50,1), col="blue", main="maxPause map3TomaxTSNs KS FDR <= 0.1",
     add=T
)

hist(sub_df_SNP$AllelicMaxPauseDist,     freq = F, ylim=c(0,0.4), las=1,
     breaks = seq(-0.5,50,1),
     col="red", density = 45, add=T)

lines(ecdf(sub_df$AllelicMaxPauseDist), col="blue")
lines(ecdf(sub_df_SNP$AllelicMaxPauseDist), col="red", pch=10)


legend("right", 
       # legend = c( paste("Control, n = ", dim(sub_df)[1], sep=""), 
       #             paste("SNP at short pause, n=", dim(sub_df_SNP)[1], sep="")),
       legend = c( paste("Allelic difference, n = ", dim(sub_df)[1], sep=""), 
                   paste("+ SNP at short pause, n=", dim(sub_df_SNP)[1], sep="")),
       title = ,
       #pch=c(15,15),
       #cex=2, 
       lty=c(0,0),
       #bty="n",
       lwd=1.5, 
       density=c(10000,25),
       angle=c(180,45),
       #angle=45,
       fill=c("blue","red"),
       border =c("black","red"),
       bty = "n"
)
legend("right", 
       # legend = c( paste("Control, n = ", dim(sub_df)[1], sep=""), 
       #             paste("SNP at short pause, n=", dim(sub_df_SNP)[1], sep="")),
       legend = c( paste("Allelic difference, n = ", dim(sub_df)[1], sep=""), 
                   paste("+ SNP at short pause, n=", dim(sub_df_SNP)[1], sep="")),
       title = ,
       #pch=c(15,15),
       #lty=c(0,0),
       #bty="n",
       lwd=1.5, 
       #density=c(10000,25),
       #angle=c(180,45),
       #angle=45,
       #fill=c("blue","dark organe","dark green")
       col=c("blue","red"),
       pch=c(19,10),
       bty = "n"
)

dev.off()

ks.test(sub_df$AllelicMaxPauseDist ,sub_df_SNP$AllelicMaxPauseDist, alternative = "less")

# indel 
sub_df_indel=df[df$map3.p.value.fdr<=0.1 & df$AllelicMaxPauseDist >0 & df$dist_Indel_maxTSN <= df$latePause,]
                #c('V1','V2','V3','earlyPause','AllelicMaxPauseDist')])
dim(sub_df_indel)
sub_df_indel$earlyPauseSites=getEarlyPauseSites(sub_df_indel)$V2
sub_df_indel=sub_df_indel[!duplicated(sub_df_indel[,c('V1','earlyPauseSites','V6','AllelicMaxPauseDist')]),]
dim(sub_df_indel)
#Figure4I
# panel56_Indel_allelicMaxPauseDist
pdf("Indel_allelicMaxPauseDist.pdf", width=7, height = 3.5, useDingbats=FALSE)
par(mar=c(6.1, 7.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
par(mgp=c(3,1,0))
par(cex.lab=2.2, cex.axis=2.2)
plot(ecdf(sub_df$AllelicMaxPauseDist), col="blue", 
     xlab="Distance between allelic maxPause",
     ylab="fraction",
     las=1,
     main="maxPause map3TomaxTSNs KS FDR <= 0.1, e and l >1bp move")

hist(sub_df$AllelicMaxPauseDist,
     freq = F, ylim=c(0,0.4), las=1,
     breaks = seq(-0.5,50,1), col="blue", main="maxPause map3TomaxTSNs KS FDR <= 0.1",
     add=T
)


hist(sub_df_indel$AllelicMaxPauseDist,     freq = F, ylim=c(0,0.4), las=1,
     breaks = seq(-0.5,50,1),
     col="dark green", density = 45, add=T)


lines(ecdf(sub_df$AllelicMaxPauseDist), cex=1, col="blue", pch=19)
lines(ecdf(sub_df_indel$AllelicMaxPauseDist), cex=1, col="dark green", pch=10)


legend("right", 
       # legend = c( paste("All FDR<=0.1, >1bp move, n = ", dim(sub_df)[1], sep=""), 
       #             paste("With Indel, n=", dim(sub_df_indel)[1], sep="")),
       legend = c( paste("Allelic difference, n = ", dim(sub_df)[1], sep=""), 
                   paste("+ Indel, n=", dim(sub_df_indel)[1], sep="")),
       title = ,
       lwd=1.5, 
       col=c("blue","dark green"),
       pch=c(19,10),
       bty = "n"
)

legend("right", 
       #legend = c( paste("All FDR<=0.1, >1bp move, n = ", dim(sub_df)[1], sep=""), 
       #            paste("Without Indel, n=", dim(sub_df_no_indel)[1], sep=""),                   paste("With Indel, n=", dim(sub_df_indel)[1], sep="")),
       legend = c( paste("Allelic difference, n = ", dim(sub_df)[1], sep=""), 
                   paste("+ Indel, n=", dim(sub_df_indel)[1], sep="")),
       
       title = ,
       #pch=c(15,15),
       #cex=1.5, 
       lty=c(0,0),
       #bty="n",
       lwd=1.5, 
       density=c(10000,25),
       angle=c(180,45),
       #angle=45,
       fill=c("blue","dark green")
       , bty = "n"
)
dev.off()

# KS test 
ks.test(sub_df$AllelicMaxPauseDist ,sub_df_SNP$AllelicMaxPauseDist, alternative = "less")
ks.test(sub_df$AllelicMaxPauseDist ,sub_df_indel$AllelicMaxPauseDist, alternative = "greater")
#ks.test(sub_df_no_indel$AllelicMaxPauseDist ,sub_df_indel$AllelicMaxPauseDist, alternative = "greater")

# count SNPs 
# use UniqRowCount dim(unique(data.frame(subdf$V2, subdf$earlyPause)))[1]
SNP_count_and_or <- function(snp_position, subDF=df){
  s1<-NULL
  s2<-NULL
  s3<-NULL
  for (i in 1:NROW(subDF)){
    SNP_distance_earlyPause = as.numeric(unlist(strsplit(as.character(subDF$SNP_distance_earlyPause[i]), ",")))
    if (length(snp_position) ==1){
      s1= c(s1, (snp_position[1]%in% SNP_distance_earlyPause))
    }
        if (length(snp_position) >=2){
      s1= c(s1, (snp_position[1]%in% SNP_distance_earlyPause))
      s2= c(s2, (snp_position[2]%in% SNP_distance_earlyPause))
    }
    if (length(snp_position) ==3){
      s3= c(s3, (snp_position[3]%in% SNP_distance_earlyPause))
    }
  }
  if (length(snp_position) == 2){
    a = s1 & s2
    o = s1 | s2
  }
  if (length(snp_position) == 3){
    a = s1 & s2 & s3
    o = s1 | s2 | s3
  }
  if (length(snp_position) == 1){
    return (list(snp_position, UniqRowCount(subDF[s1,])))
  }
  if (length(snp_position) == 2){
    return (list(UniqRowCount(subDF[a,]),UniqRowCount(subDF[o,]),
                 snp_position, UniqRowCount(subDF[s1,]), UniqRowCount(subDF[s2,])))
  }
    if (length(snp_position) >= 3){
  return (list(UniqRowCount(subDF[a,]),UniqRowCount(subDF[o,]),
               snp_position, UniqRowCount(subDF[s1,]), UniqRowCount(subDF[s2,]), UniqRowCount(subDF[s3,])))
}
}

# Figure4F
df_target=df[df$map3.p.value.fdr<=0.1&  # KS test 
               df$AllelicMaxPauseDist >0,] # at least 1bp apart]
dim(df_target)
UniqRowCount(df_target)


# pause distribution not significantly different
df_target=df[df$map3.p.value.fdr> 0.9 & # KS test 
               df$AllelicMaxPauseDist ==0,] # at least 1bp apart]
sum(df$map3.p.value.fdr> 0.9)
dim(df_target)
UniqRowCount(df_target)

SNP_count_and_or(c(0,-2),df_target) 
SNP_count_and_or(c(0,-3),df_target)
SNP_count_and_or(c(-2, -3),df_target)
SNP_count_and_or(c(1,0),df_target)
SNP_count_and_or(c(1:10),df_target)
SNP_count_and_or(c(1,-3),df_target)
SNP_count_and_or(c(1, -2),df_target)
SNP_count_and_or(c(1, -2, -3),df_target)

# with indel
SNP_count_and_or(0, subDF = df_target[df_target$dist_Indel_maxTSN <= df_target$latePause,])
SNP_count_and_or(c(-2,-3), subDF = df_target[df_target$dist_Indel_maxTSN <= df_target$latePause,])
SNP_count_and_or(-3, subDF = df_target[df_target$dist_Indel_maxTSN <= df_target$latePause,])
SNP_count_and_or(-2, subDF = df_target[df_target$dist_Indel_maxTSN <= df_target$latePause,])
SNP_count_and_or(1, subDF = df_target[df_target$dist_Indel_maxTSN <= df_target$latePause,])

# without indel
SNP_count_and_or(0, subDF = df_target[df_target$dist_Indel_maxTSN > df_target$latePause,])
SNP_count_and_or(-2, subDF = df_target[df_target$dist_Indel_maxTSN > df_target$latePause,])
SNP_count_and_or(-3, subDF = df_target[df_target$dist_Indel_maxTSN > df_target$latePause,])
SNP_count_and_or(c(-2,-3), subDF = df_target[df_target$dist_Indel_maxTSN > df_target$latePause,])
SNP_count_and_or(1, subDF = df_target[df_target$dist_Indel_maxTSN > df_target$latePause,])

# examine if there is SNPs at early pause in regions with fdr>0.9 # pause distribution not significantly different
sub_df_ns=df[df$map3.p.value.fdr > 0.9 & df$AllelicMaxPauseDist ==0,]
sub_df_SNP_ns=df[df$map3.p.value.fdr > 0.9 & df$AllelicMaxPauseDist ==0  & df$SNP1_distance_earlyPause==0 & !is.na(df$SNP1_distance_earlyPause),]
#View(sub_df_SNP_ns)
UniqRowCount(sub_df_ns)  #1396
UniqRowCount(sub_df_SNP_ns) # 29
# examine if the SNPs at C were enriched in pause region with difference (KS test, fdr<=0.1)
fisher.test(data.frame(
  UniqRowCount(sub_df_SNP),UniqRowCount(sub_df),
    c(UniqRowCount(sub_df_SNP_ns),UniqRowCount(sub_df_ns))))
#p-value < 0.00000000000000022
#####HERE#####

### What kind of SNPs?

df$earlyPause_parent = "M"
df$earlyPause_parent[df$earlyPause==df$pat_maxPause_map3] = "P"
df$earlyPause_parent[df$earlyPause==df$mat_maxPause_map3] = "M"
sum(df$earlyPause_parent=="M"); sum(df$earlyPause_parent=="P")
show.window=0
# early pause, background 
bed0 <- df[,1:6]
bed0$earlyPause_parent = df$earlyPause_parent
bed0$V4 = df$earlyPause_parent
show.window=0
for (i in 1:NROW(bed0)){
  if(bed0[i,6]=="-") {
    bed0[i,3] <- df[i,3] - df$earlyPause[i] + show.window
    bed0[i,2] <- df[i,2] - df$earlyPause[i] - show.window
  } else {
    bed0[i,2] <- df[i,2] + df$earlyPause[i] - show.window
    bed0[i,3] <- df[i,3] + df$earlyPause[i] + show.window
  }
}
bed0 = unique(bed0[df$map3.p.value.fdr>0.9 & (df$mat_maxPause_map3 == df$pat_maxPause_map3) ,])
dim(bed0) #1396
write.table(bed0, file="Tissues3_EarlyPause_BG.bed", quote = F, sep="\t", row.names = F, col.names = F)


bed7 <- df[,1:6]
bed7$V4 = df$earlyPause_parent
bed7$earlyPause_parent = df$earlyPause_parent
for (i in 1:NROW(bed7)){
  if(bed7[i,6]=="-") {
    bed7[i,3] <- df[i,3] - df$earlyPause[i] + show.window
    bed7[i,2] <- df[i,2] - df$earlyPause[i] - show.window
  } else {
    bed7[i,2] <- df[i,2] + df$earlyPause[i] - show.window
    bed7[i,3] <- df[i,3] + df$earlyPause[i] + show.window
  }
}
# early and late pause with at least 1 bp apart AND map3 KS test fdr<=0.1
#early
bed7 = bed7[df$earlyPause != df$latePause & df$map3.p.value.fdr<=0.1,]
bed7 = unique(bed7)
dim(bed7)
dim(df)
#Figure 4G, left
write.table(bed7, file="Tissues3_EarlyPause_1bpapart_KSfdr0.1.bed", quote = F, sep="\t", row.names = F, col.names = F)



bed9 <- df[,1:6]
bed9$earlyPause_parent = df$earlyPause_parent
bed9$V4 = df$earlyPause_parent
for (i in 1:NROW(bed9)){
  if(bed9[i,6]=="-") {
    bed9[i,3] <- df[i,3] - df$latePause[i] 
    bed9[i,2] <- df[i,2] - df$latePause[i] 
  } else {
    bed9[i,2] <- df[i,2] + df$latePause[i]  
    bed9[i,3] <- df[i,3] + df$latePause[i]  
  }
}
# early and late pause with at least 1 bp apart AND map3 KS test fdr<=0.1
# late
bed9 = unique(bed9[df$earlyPause != df$latePause & df$map3.p.value.fdr<=0.1,])
dim(bed9)
dim(df)
#Figure 4G, right
write.table(bed9, file="Tissues3_LatePause_1bpapart_KSfdr0.1.bed", quote = F, sep="\t", row.names = F, col.names = F)


# maxPause from all (mat, pat, ide) reads, combine all organs, duplicates removed 
bed8 <- df[,1:6]
for (i in 1:NROW(bed8)){
  if(bed8[i,6]=="-") {
    bed8[i,3] <- df[i,3] - df$maxPause_map3AllReads[i]
    bed8[i,2] <- df[i,2] - df$maxPause_map3AllReads[i]
  } else {
    bed8[i,2] <- df[i,2] + df$maxPause_map3AllReads[i]
    bed8[i,3] <- df[i,3] + df$maxPause_map3AllReads[i]
  }
}

dim(bed8)
bed8 = bed8[!duplicated(bed8$V2),]
#Figure 4D bottom
dim(bed8)
write.table(bed8, file="combine_maxPause_noduplicate.bed", quote = F, sep="\t", row.names = F, col.names = F)



library("TmCalculator")
library(seqLogo)
#define function that divides the frequency by the row sum i.e. proportions
proportion <- function(x){
  rs <- sum(x);
  return(x / rs);
}

seq_upperCase <- function(seq){
  seq[seq=="a"]<- "A"
  seq[seq=="t"]<- "T"
  seq[seq=="c"]<- "C"
  seq[seq=="g"]<- "G"
  
  a <- NULL
  t <- NULL
  c <- NULL
  g <- NULL
  for (i in 1:NCOL(seq)){
    a <- c(a, sum(seq[,i]=="A"))
    t <- c(t, sum(seq[,i]=="T"))
    c <- c(c, sum(seq[,i]=="C"))
    g <- c(g, sum(seq[,i]=="G"))
  }
  return (data.frame(a,c,g,t))
}

SeqLogo <- function(seq, output, range=NULL) {
  #seq=m$HighAlleleSeq
  seq<- data.frame(do.call(rbind, strsplit(as.character(seq), "")))
  df <- seq_upperCase(seq)
  #create position weight matrix
  pwm <- apply(df, 1, proportion)
  if (!is.null(range)){
    p = makePWM((pwm[,range]))    
  }else{
    p = makePWM((pwm))  
  }
  
  #p <- makePWM(pwm)
  # slotNames(p)
  # p@consensus
  # p@ic
  # p@width
  # p@alphabet
  pdf(output)
  seqLogo(p)
  dev.off()
  return (pwm)
}

# Panel58_short_pause_site_seqlogo
seq_a=read.table("Tissues3_EarlyPause_1bpapart_KSfdr0.1_+-10_Early_LateAlleleSeq.bed")
dim(seq_a)
#Figure4G left, n=270
Tissues3_EarlyPause_1bpapart_KSfdr0.1_early=SeqLogo(seq_a$V8, "Tissues3_EarlyPause_1bpapart_KSfdr0.1_early.pdf")
Tissues3_EarlyPause_1bpapart_KSfdr0.1_late=SeqLogo(seq_a$V9, "Tissues3_EarlyPause_1bpapart_KSfdr0.1_late.pdf")

seq_l=read.table("Tissues3_LatePause_1bpapart_KSfdr0.1_+-10_Early_LateAlleleSeq.bed")
dim(seq_l)
#Figure4G right, n=278
Tissues3_LatePause_1bpapart_KSfdr0.1_early=SeqLogo(seq_l$V8, "Tissues3_LatePause_1bpapart_KSfdr0.1_early.pdf")
Tissues3_LatePause_1bpapart_KSfdr0.1_late=SeqLogo(seq_l$V9, "Tissues3_LatePause_1bpapart_KSfdr0.1_late.pdf")

seq=read.table("combine_maxPause_noduplicate_+-30_mm10_Seq.bed")
dim(seq)
#Figure4D bottom
#n=3456
combine_maxPause_noduplicate= SeqLogo(seq$V7, "combine_maxPause_noduplicate_+-30_mm10_Se.pdf")


acgt_col=c("dark green", "blue", "orange" , "red")
acgt=c("A","C","G","T")

# combine maxPause 
C <- function(a){
  return (sum(a=="C"|a=="c")/length(a)*100)
}
G <- function(a){
  return (sum(a=="G"|a=="g")/length(a)*100)
}
par(mfrow=c(1,3))
bin=10
d=30
#for (bin in c(10,12)){
seq_df=read.table("combine_maxPause_noduplicate_+-30_mm10_Seq.bed")
dim(seq_df)
range1=(d-bin*2+1):(d-bin)
range2=(d-bin+1):d
range3=(d+2):(d+bin)

for (i in 1:NROW(seq)){
  a=s2c(as.character(seq$V7[i]))
  a1=a[range1]
  a2=a[range2]
  a3=a[range3]
  seq_df$GC_range1[i]= GC(a1)
  seq_df$GC_range2[i]= GC(a2)
  seq_df$GC_range3[i]= GC(a3)
  seq_df$G_range1[i]= G(a1)
  seq_df$G_range2[i]= G(a2)
  seq_df$G_range3[i]= G(a3)
  seq_df$C_range1[i]= C(a1)
  seq_df$C_range2[i]= C(a2)
  seq_df$C_range3[i]= C(a3)
  #sum(a=="C"|a=="c")/length(a)*100
}

library("vioplot")
# panel57_GC_content_around_maxPause_noduplicate
# Figure4E, n=3456
dim(seq_df)
pdf("GC_content_around_maxPause_noduplicate.pdf", width=7, height = 3.5, useDingbats=FALSE)
par(mfrow=c(1,3))
par(mar=c(6.1, 7.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
par(mgp=c(3,1,0))
par(cex.lab=2.2, cex.axis=2.2)
vioplot(seq_df$G_range1, seq_df$G_range2, seq_df$G_range3, 
        main= paste("bin size = ", bin, sep=""),
        ylab="G content", ylim=c(0,100), las=1, frame.plot = F)
#abline(h=40, col="yellow")
vioplot(seq_df$GC_range1, seq_df$GC_range2, seq_df$GC_range3, 
        main= paste("bin size = ", bin, sep=""),
        ylab="GC content", ylim=c(0,100), las=1, frame.plot = F)
#abline(h=60, col="green")
vioplot(seq_df$C_range1, seq_df$C_range2, seq_df$C_range3, 
        main= paste("bin size = ", bin, sep=""),
        ylab="C content", ylim=c(0,100), las=1, frame.plot = F)
#abline(h=20, col="yellow")
dev.off()

# Wilcoxon Rank Sum and Signed Rank Tests
wilcox.test(seq_df$G_range1, seq_df$G_range2, paired = TRUE)
wilcox.test(seq_df$G_range2, seq_df$G_range3, paired = TRUE)
wilcox.test(seq_df$G_range1, seq_df$G_range3, paired = TRUE)

wilcox.test(seq_df$GC_range1, seq_df$GC_range2, paired = TRUE)
wilcox.test(seq_df$GC_range2, seq_df$GC_range3, paired = TRUE)
wilcox.test(seq_df$GC_range1, seq_df$GC_range3, paired = TRUE)

wilcox.test(seq_df$C_range1, seq_df$C_range2, paired = TRUE)
wilcox.test(seq_df$C_range2, seq_df$C_range3, paired = TRUE)
wilcox.test(seq_df$C_range1, seq_df$C_range3, paired = TRUE)

# early background
seq_bg=read.table("Tissues3_EarlyPause_BG_+-10_Early_LateAlleleSeq.bed")
dim(seq_bg)
Tissues3_EarlyPause_BG_early=SeqLogo(seq_bg$V8, "Tissues3_EarlyPause_BG_early.pdf")
Tissues3_EarlyPause_BG_late=SeqLogo(seq_bg$V9, "Tissues3_EarlyPause_BG_late.pdf")

#Sup Figure 4B
pdf(paste(organ,"_pause_deltaATCG.pdf" ,sep=""), width =7, height = 7,useDingbats=FALSE)
par(mfcol=c(3,1))
par(mar=c(6.1, 7.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
par(mgp=c(3,1,0))
par(cex.lab=2.2, cex.axis=2.2)
#par(cex=1.5)
pch_u=1
w=10


organ="3 Tissues" ; target=Tissues3_EarlyPause_BG_early - Tissues3_EarlyPause_BG_late
plot(-w:w, target[1,] , col = acgt_col[1], type="o",
     ylim=c(-0.15,0.15), 
     pch=pch_u,
     ylab="Short Allele - Long Allele",
     xlab="Distance to short allelic maxPause",
     main=paste(organ,"BG, n=",dim(seq_bg)[1], sep=" "),
     las=1, frame=FALSE
)
abline(h=0, col="gray")
for (i in 4:1){
  points(-w:w,target[i,], col = acgt_col[i], type="o", pch=pch_u)
}
legend("topleft", legend=acgt,
       col=acgt_col, bty = "n", lty=1, pch=pch_u, cex=1.5)

#target=test_early - test_late
organ="3 Tissues" ; target=Tissues3_EarlyPause_1bpapart_KSfdr0.1_early - Tissues3_EarlyPause_1bpapart_KSfdr0.1_late
plot(-w:w, target[1,] , col = acgt_col[1], type="o",
     ylim=c(-0.15,0.15), 
     pch=pch_u,
     ylab="Short Allele - Long Allele",
     xlab="Distance to short allelic maxPause",
     main=paste(organ,"E L at least 1bp apart, KS test fdr<=0.1, n=",dim(seq_a)[1], sep=" "),
     las=1, frame=FALSE
)
abline(h=0, col="gray")
for (i in 4:1){
  points(-w:w,target[i,], col = acgt_col[i], type="o", pch=pch_u)
}
#legend("topleft", legend=acgt,col=acgt_col, bty = "n", lty=1, pch=pch_u)

organ="3 Tissues" ; target=Tissues3_LatePause_1bpapart_KSfdr0.1_early - Tissues3_LatePause_1bpapart_KSfdr0.1_late
plot(-w:w, target[1,] , col = acgt_col[1], type="o",
     ylim=c(-0.15,0.15), 
     pch=pch_u,
     ylab="Short Allele - Long Allele",
     xlab="Distance to long allelic maxPause",
     main=paste(organ,"E L at least 1bp apart, KS test fdr<=0.1, n=", dim(seq_l)[1], sep=" "),
     las=1, frame=FALSE
)
abline(h=0, col="gray")
for (i in 4:1){
  points(-w:w,target[i,], col = acgt_col[i], type="o", pch=pch_u)
}
#legend("topleft", legend=acgt,
#       col=acgt_col, bty = "n", lty=1, pch=pch_u)

dev.off()

