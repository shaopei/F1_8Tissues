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


# get read length distribution
for (Tissue in c("HT","SK", "KD")){
  pdf(file = paste(Tissue, "_read_length_distribution.pdf"))
  ide_name=paste(Tissue,"_mat_identical_read_length.txt",sep="")
  mat_name=paste(Tissue,"_mat_specific_read_length.txt",sep="")
  pat_name=paste(Tissue,"_pat_specific_read_length.txt",sep="")
  
  ide<-read.table(ide_name)
  mat<-read.table(mat_name)
  pat<-read.table(pat_name)
  
  par(mfrow=c(3,1))
  xlim_top=101
  hist(mat$V1[mat$V1<=100], 
       xlim = c(0, xlim_top), las=1, 
       breaks = seq(0.5, xlim_top,1),
       main=Tissue)
  hist(pat$V1[pat$V1<=100], 
       breaks = seq(0.5, xlim_top,1),
       xlim = c(0, xlim_top), las=1)
  hist(ide$V1[ide$V1<=100], 
       breaks = seq(0.5, xlim_top,1),
       xlim = c(0, xlim_top), las=1)
  dev.off()
}


# pause center analysis, using dREG sites with KS test, map3 position, fdr<=0.1
combine_pause0.1 <- NULL
#t="KD"
for (t in c("HT","SK","KD")){
  show.window=100
  pause_window_0.1 <- read.table(paste(file_dir,t,"_dREG_5mat5pat_uniq_pValue_fdr0.1.bed", sep =""), header = F)
  pause_window_0.1 <- pause_window_0.1[pause_window_0.1$V1 != 'chrX',]
  end=".rpm.bw"; times=10
  #end=".bw"; times=1
  # allelic reads (map3)
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

old_AT = AT
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
hist(AT$pause.dist[AT$TSN.dist==0]
     , breaks=seq(-0.5,35,1)
)
hist(abs(AT$TSN.dist[AT$pause.dist==0])# &AT$TSN.dist !=0 ]
     , breaks=seq(-0.5,35,1)
     #, ylim=c(0,20)
)
sum(AT$pause.dist==0 &AT$TSN.dist !=0)

View(AT[AT$pause.dist==0 &AT$TSN.dist !=0, ])

sum(AT$TSN.dist==0)
sum(AT$TSN.dist!=0)
#hist(AT$pause.dist[AT$TSN.dist != 0] - AT$TSN.dist[AT$TSN.dist != 0]
#     , breaks=seq(-35.5,35,1)
#)


subAT=AT[abs(AT$pause.dist - AT$TSN.dist) > 10,]
subAT=subAT[subAT$TSN.dist==0,]
subAT$deltaTsnPauseDist = subAT$pause.dist - subAT$TSN.dist

library("ggplot2")

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


for (i in 1:dim(df)[1]){
  m=as.numeric(unlist(strsplit(as.character(df$V7[i]), ",")))
  p=as.numeric(unlist(strsplit(as.character(df$V8[i]), ",")))
  df$RL.p.value[i] = ks.test(m,p) $ p.value
}
for (i in 1:dim(df)[1]){
  m=as.numeric(unlist(strsplit(as.character(df$V9[i]), ",")))
  p=as.numeric(unlist(strsplit(as.character(df$V10[i]), ",")))
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

UniqRowCount <- function(subdf){
 dim(unique(data.frame(subdf$V2, subdf$earlyPause)))[1]
}


df$target=df$V15 <= df$maxPauseSite_map3 & (df$map3.p.value.fdr <=0.1 ) & (df$mat_maxPause_map3 != df$pat_maxPause_map3)
dim(df)
# use only fdr<=0.1 sites where early pause != late pause
sum((df$map3.p.value.fdr <=0.1 ) & (df$mat_maxPause_map3 != df$pat_maxPause_map3))
a=UniqRowCount(df[(df$map3.p.value.fdr <=0.1 ) & (df$mat_maxPause_map3 != df$pat_maxPause_map3),])
b=UniqRowCount(df[df$V15 <= df$maxPauseSite_map3 & (df$map3.p.value.fdr <=0.1 ) & (df$mat_maxPause_map3 != df$pat_maxPause_map3),])
c=UniqRowCount(df[(df$map3.p.value.fdr >0.9 ) & (df$mat_maxPause_map3 == df$pat_maxPause_map3),])
d=UniqRowCount(df[df$V15 <= df$maxPauseSite_map3 &(df$map3.p.value.fdr >0.9 ) & (df$mat_maxPause_map3 == df$pat_maxPause_map3),])
a;b;c;d
fisher.test(data.frame(c(b, a), c(d, c)))
subdf=(df[(df$map3.p.value.fdr <=0.1 ) & (df$mat_maxPause_map3 != df$pat_maxPause_map3),c(1,2,3,14,5,6)])
dim(unique(subdf))

#sum(df$target)
#names(df)
#subdf = (df[(df$map3.p.value.fdr <=0.1 ) & (df$mat_maxPause_map3 != df$pat_maxPause_map3), c(1:6, 25,26, 19,20,16)])
#sum(df$map3.p.value.fdr<=0.1)
#sum(df$map3.p.value.fdr>0.1)

#sum(df$V15 > df$maxPauseSite_map3 & (df$map3.p.value.fdr <=0.1 ))
#sum(df$V15 > df$maxPauseSite_map3 & (df$map3.p.value.fdr > 0.1 ))

# plot( (df$mat_AvePause_map3- df$pat_AvePause_map3)[df$target], 
#       (df$mat_Idel_length - df$pat_Idel_length)[df$target],
#       xlim=c(-1*lim,lim), ylim=c(-1*lim,lim),
#       xlab="B6 - CAST pause position",
#       ylab="B6 - CAST indel length",
#       frame=F, 
#       pch=19, col=rgb(0,0,0,alpha = 0.25))
pch_u=19
plot( (df$mat_AvePause_map3- df$pat_AvePause_map3)[df$target& df$Tissue=="HT"],
      (df$mat_Idel_length - df$pat_Idel_length)[df$target & df$Tissue=="HT"],
      xlim=c(-1*lim,lim),
      ylim=c(-1*lim,lim),
      xlab="B6 - CAST pause position",
      ylab="B6 - CAST indel length",
      pch=pch_u, frame=F,
      col="red")
# points( (df$mat_AvePause_map3- df$pat_AvePause_map3)[df$target& df$Tissue=="HT"], 
#         (df$mat_Idel_length - df$pat_Idel_length)[df$target & df$Tissue=="HT"],
#         pch=pch_u, 
#         col="red")
points( (df$mat_AvePause_map3- df$pat_AvePause_map3)[df$target& df$Tissue=="SK"], 
        (df$mat_Idel_length - df$pat_Idel_length)[df$target & df$Tissue=="SK"],
        pch=pch_u, 
        col="dark orange")
points( (df$mat_AvePause_map3- df$pat_AvePause_map3)[df$target& df$Tissue=="KD"], 
        (df$mat_Idel_length - df$pat_Idel_length)[df$target & df$Tissue=="KD"],
        pch=pch_u, col="blue")

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
       col=c("red", "dark orange", "blue")
       , bty = "n"
)



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
  bed6 = unique(bed6)
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
tempFunc <-function(bp.apart=1, upto=100){
  step=1
  
  ## SNPs around early pause
  t1=metaplot.SNPsLocation.aroundEarlyPause(df[df$map3.p.value.fdr<=0.1 
                                            & abs(df$mat_maxPause_map3 - df$pat_maxPause_map3) >= bp.apart 
                                            & abs(df$mat_maxPause_map3 - df$pat_maxPause_map3) <= upto
                                            & df$dist_Indel_maxTSN > df$latePause,], 
                                         name=paste("HSK, step=", step,", early and late at least ", bp.apart," bp apart, upto " ,upto, sep=""), col="red", step=step)
  abline(v=0, col="gray")

  t2=metaplot.SNPsLocation.aroundEarlyPause(df[df$map3.p.value.fdr>0.9  
                                               & (df$mat_maxPause_map3 == df$pat_maxPause_map3) 
                                            & df$dist_Indel_maxTSN > df$latePause,], 
                                         col="blue", step=step, add=TRUE)

  legend("topright", legend=c(paste("Early Pause, fdr <= 0.1, n= ", t1[[3]], sep=""),
                              paste("Early Pause, fdr >  0.9, n= ", t2[[3]], sep="")),
         col=c("red", "blue"), 
         bty = "n", lty=1, pch=19)

}

# fdr<=0.1, early != late pause,  and with indel
a=UniqRowCount(df[(df$map3.p.value.fdr <=0.1 ) & (df$mat_maxPause_map3 != df$pat_maxPause_map3),])
b=UniqRowCount(df[df$dist_Indel_maxTSN <= df$latePause & (df$map3.p.value.fdr <=0.1 ) & (df$mat_maxPause_map3 != df$pat_maxPause_map3),])
c=UniqRowCount(df[(df$map3.p.value.fdr >0.9 ) & (df$mat_maxPause_map3 == df$pat_maxPause_map3),])
d=UniqRowCount(df[df$dist_Indel_maxTSN <= df$latePause &(df$map3.p.value.fdr >0.9 ) & (df$mat_maxPause_map3 == df$pat_maxPause_map3),])
a;b;c;d
fisher.test(data.frame(c(b, a), c(d, c)))

if(0){
sum(df$map3.p.value.fdr<=0.1)
sum(df$earlyPause==df$latePause)


sum(df$map3.p.value.fdr<=0.1 & df$dist_Indel_maxTSN <= df$latePause)

sum(df$map3.p.value.fdr<=0.1 & abs(df$mat_maxPause_map3 - df$pat_maxPause_map3) ==0
    & df$dist_Indel_maxTSN <= df$latePause)
sum(df$map3.p.value.fdr<=0.1 & abs(df$mat_maxPause_map3 - df$pat_maxPause_map3) >= bp.apart & abs(df$mat_maxPause_map3 - df$pat_maxPause_map3) <= upto)
sum(df$map3.p.value.fdr<=0.1 & abs(df$mat_maxPause_map3 - df$pat_maxPause_map3) >= bp.apart & abs(df$mat_maxPause_map3 - df$pat_maxPause_map3) <= upto
    & df$dist_Indel_maxTSN <= df$latePause)


subdf=df[df$map3.p.value.fdr<=0.1 & abs(df$mat_maxPause_map3 - df$pat_maxPause_map3) >= bp.apart & abs(df$mat_maxPause_map3 - df$pat_maxPause_map3) <= upto,]
dim(subdf)
sum(subdf$dist_Indel_maxTSN <= subdf$latePause)
}

##
## SNPs around early pause, use.sum = TRUE, for test
step=1; bp.apart=1; upto=100
e1=metaplot.SNPsLocation.aroundEarlyPause(df[df$map3.p.value.fdr<=0.1 
                                             & abs(df$mat_maxPause_map3 - df$pat_maxPause_map3) >= bp.apart 
                                             & abs(df$mat_maxPause_map3 - df$pat_maxPause_map3) <= upto
                                             & df$dist_Indel_maxTSN > df$latePause,], 
                                          name=paste("HSK, step=", step,", early and late at least ", bp.apart," bp apart, upto " ,upto, sep=""), col="red", step=step,
                                          use.sum = TRUE)
abline(v=0, col="gray")
# abline(v=-5, col="gray")
# abline(v=-10, col="gray")

e9=metaplot.SNPsLocation.aroundEarlyPause(df[df$map3.p.value.fdr>0.9
                                             & (df$mat_maxPause_map3 == df$pat_maxPause_map3) 
                                             & df$dist_Indel_maxTSN > df$latePause,], col="blue", step=step, use.sum = TRUE)
# legend("topright", legend=c(paste("fdr <= 0.1, n= ", sum(df$map3.p.value.fdr<=0.1 & abs(df$mat_maxPause_map3 - df$pat_maxPause_map3) >0), sep=""), 
#                             paste("fdr >  0.9, n= ", sum(df$map3.p.value.fdr>0.9), sep="")),
#        col=c("red", "blue"), bty = "n", lty=1, pch=19)
# 

# add a test
#pdf(paste(organ,"_FisherExactTest_CountOfTSSwithAT2GCSNPs_vs_AllTSS.pdf", sep=""))
par(mfrow=c(3,1))
tempFunc(1,100) # use mean
fdr_cutoff=0.1
if(0){
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


# effect of SNPs at the difference between maxPauses of early and late pause
# exclude allelic maxPause difference 0
df$AllelicMaxPauseDist = abs(df$mat_maxPause_map3 - df$pat_maxPause_map3)

sub_df=df[df$map3.p.value.fdr<=0.1 & df$AllelicMaxPauseDist >0,]
sub_df_SNP=df[df$map3.p.value.fdr<=0.1 & df$AllelicMaxPauseDist >0 & 
                df$SNP1_distance_earlyPause==0 & !is.na(df$SNP1_distance_earlyPause),]
remove_duplicates_V2_earlyPause <- function(subdf){
  return(subdf[!duplicated(subdf[,c('V2','earlyPause')]),])
}

set_remove_duplicates_V2_earlyPause =FALSE
if (set_remove_duplicates_V2_earlyPause){
  sub_df=remove_duplicates_V2_earlyPause(sub_df)
  sub_df_SNP=remove_duplicates_V2_earlyPause(sub_df_SNP)
}else{

  sub_df=unique(df[df$map3.p.value.fdr<=0.1 & df$AllelicMaxPauseDist >0, 
                   c('V1','V2','V3','earlyPause','AllelicMaxPauseDist')])
  sub_df_SNP=unique(df[df$map3.p.value.fdr<=0.1 & df$AllelicMaxPauseDist >0 
                     & df$SNP1_distance_earlyPause==0 & !is.na(df$SNP1_distance_earlyPause),
                     c('V1','V2','V3','earlyPause','AllelicMaxPauseDist')])
}

#par(mfrow=c(2,1))
plot(ecdf(sub_df$AllelicMaxPauseDist), col="blue", 
     xlab="AllelicMaxPauseDist",
     ylab="density",
     las=1,
     main="maxPause map3TomaxTSNs KS FDR <= 0.1, e and l >1bp move")
hist(sub_df$AllelicMaxPauseDist,
     freq = F, ylim=c(0,0.4), las=1,
     breaks = seq(-0.5,50,1), col="blue", main="maxPause map3TomaxTSNs KS FDR <= 0.1",
     add=T
)

hist(sub_df_SNP$AllelicMaxPauseDist,     freq = F, ylim=c(0,0.4), las=1,
     breaks = seq(-0.5,50,1),
     col="red", density = 45, add=T)

lines(ecdf(sub_df$AllelicMaxPauseDist), col="blue")
lines(ecdf(sub_df_SNP$AllelicMaxPauseDist), col="red")


legend("right", 
       legend = c( paste("All FDR<=0.1, >1bp move, n = ", dim(sub_df)[1], sep=""), 
                   paste("SNP at early pause, n=", dim(sub_df_SNP)[1], sep="")),
       title = ,
       #pch=c(15,15),
       #cex=2, 
       lty=c(0,0),
       #bty="n",
       lwd=1.5, 
       density=c(10000,25),
       angle=c(180,45),
       #angle=45,
       fill=c("blue","red")
       , bty = "n"
)

ks.test(sub_df$AllelicMaxPauseDist ,sub_df_SNP$AllelicMaxPauseDist, alternative = "less")

# indel 
sub_df_indel=unique(df[df$map3.p.value.fdr<=0.1 & df$AllelicMaxPauseDist >0 & df$dist_Indel_maxTSN <= df$latePause,
                c('V1','V2','V3','earlyPause','AllelicMaxPauseDist')])
#sub_df_no_indel=df[df$map3.p.value.fdr<=0.1 & df$AllelicMaxPauseDist >0 & df$dist_Indel_maxTSN > df$latePause,]
#sub_df_SNP_indel=unique(df[df$map3.p.value.fdr<=0.1 & df$AllelicMaxPauseDist >0 & df$SNP1_distance_earlyPause==0 & !is.na(df$SNP1_distance_earlyPause)& df$dist_Indel_maxTSN <= df$latePause,]
dim(sub_df_indel)
#dim(sub_df_SNP_indel)

plot(ecdf(sub_df$AllelicMaxPauseDist), col="blue", 
     xlab="AllelicMaxPauseDist",
     ylab="density",
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


lines(ecdf(sub_df$AllelicMaxPauseDist), col="blue")
lines(ecdf(sub_df_indel$AllelicMaxPauseDist), col="dark green")


legend("right", 
       legend = c( paste("All FDR<=0.1, >1bp move, n = ", dim(sub_df)[1], sep=""), 
                   #paste("Without Indel, n=", dim(sub_df_no_indel)[1], sep=""),
                   paste("With Indel, n=", dim(sub_df_indel)[1], sep="")),
       title = ,
       #pch=c(15,15),
       #cex=2, 
       #lty=c(0,0),
       #bty="n",
       lwd=1.5, 
       #density=c(10000,25),
       #angle=c(180,45),
       #angle=45,
       #fill=c("blue","dark organe","dark green")
       col=c("blue","dark green"),
       pch=16,
       bty = "n"
)


# KS test 
ks.test(sub_df$AllelicMaxPauseDist ,sub_df_SNP$AllelicMaxPauseDist, alternative = "less")
ks.test(sub_df$AllelicMaxPauseDist ,sub_df_indel$AllelicMaxPauseDist, alternative = "greater")
#ks.test(sub_df_no_indel$AllelicMaxPauseDist ,sub_df_indel$AllelicMaxPauseDist, alternative = "greater")

# count SNPs 
# use UniqRowCount dim(unique(data.frame(subdf$V2, subdf$earlyPause)))[1]
SNP_count_and_or <- function(snp_position, subDF=df){
  a = (subDF[((subDF$SNP1_distance_earlyPause %in% snp_position & !is.na(subDF$SNP1_distance_earlyPause)) 
              #&(subDF$SNP3_distance_earlyPause %in% snp_position & !is.na(subDF$SNP3_distance_earlyPause)) 
              & (subDF$SNP2_distance_earlyPause %in% snp_position & !is.na(subDF$SNP2_distance_earlyPause))),])
  o = 
    (subDF[(subDF$SNP1_distance_earlyPause %in% snp_position & !is.na(subDF$SNP1_distance_earlyPause)) 
           | (subDF$SNP2_distance_earlyPause %in% snp_position & !is.na(subDF$SNP2_distance_earlyPause))
           | (subDF$SNP3_distance_earlyPause %in% snp_position & !is.na(subDF$SNP3_distance_earlyPause)),])
  count <-NULL
  for (s in snp_position){
    count = c(count, 
              UniqRowCount(subDF[((subDF$SNP1_distance_earlyPause==s & !is.na(subDF$SNP1_distance_earlyPause)) 
                                  |(subDF$SNP3_distance_earlyPause==s & !is.na(subDF$SNP3_distance_earlyPause)) 
                                  | (subDF$SNP2_distance_earlyPause==s & !is.na(subDF$SNP2_distance_earlyPause))),]))
  }
  
  return (list(UniqRowCount(a),UniqRowCount(o),snp_position, count))
  # return and, or, snp position, count of sites with snp in that position
}

df_target=df[df$map3.p.value.fdr<=0.1&  # KS test 
               df$AllelicMaxPauseDist >0,] # at least 1bp apart]

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
SNP_count_and_or(c(0,-2),df_target)
SNP_count_and_or(c(0,-3),df_target)
SNP_count_and_or(c(1:10),df_target)
SNP_count_and_or(c(1,-3),df_target)
SNP_count_and_or(c(-2, 1),df_target)

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
SNP_count_and_or(-3, subDF = df_target[df_target$dist_Indel_maxTSN > df_target$latePause,])
SNP_count_and_or(1, subDF = df_target[df_target$dist_Indel_maxTSN > df_target$latePause,])

#####HERE#####




sum(df$map3.p.value.fdr<=0.1)
# fdr<=0.1 and with indel
sum(df$map3.p.value.fdr<=0.1 & df$dist_Indel_maxTSN <= df$latePause)

sum(sub_df$map3.p.value.fdr<=0.1& sub_df$dist_Indel_maxTSN <= sub_df$latePause)
sum(sub_df_SNP$map3.p.value.fdr<=0.1& sub_df_SNP$dist_Indel_maxTSN <= sub_df_SNP$latePause)
sum(sub_df_SNP$map3.p.value.fdr<=0.1& sub_df_SNP$dist_Indel_maxTSN > sub_df_SNP$latePause)


# examine if there is SNPs at C of early pause in regions with fdr>0.9
sub_df_ns=df[df$map3.p.value.fdr > 0.9 & df$AllelicMaxPauseDist ==0,]
sub_df_SNP_ns=df[df$map3.p.value.fdr > 0.9 & df$AllelicMaxPauseDist ==0  & df$SNP1_distance_earlyPause==0 & !is.na(df$SNP1_distance_earlyPause),]
#View(sub_df_SNP_ns)
dim(sub_df_ns)  #1657
dim(sub_df_SNP_ns) # 31
# examine if the SNPs at C were enriched in pause region with difference (KS test, fdr<=0.1)
fisher.test(data.frame(c(45,308),c(31,1657)))



# SNPs at -3 G
# exclude Allelic Max Pause Distance  0
df$AllelicMaxPauseDist = abs(df$mat_maxPause_map3 - df$pat_maxPause_map3)
sub_df=df[df$map3.p.value.fdr<=0.1 & df$AllelicMaxPauseDist >0,]
sub_df_SNP_3G=df[df$map3.p.value.fdr<=0.1 & df$AllelicMaxPauseDist >0 &((df$SNP1_distance_earlyPause==-3 & !is.na(df$SNP1_distance_earlyPause))|(df$SNP2_distance_earlyPause==-3 & !is.na(df$SNP2_distance_earlyPause))),]
par(mfrow=c(2,1))
plot(ecdf(sub_df$AllelicMaxPauseDist), col="blue", 
     xlab="AllelicMaxPauseDist",
     ylab="density",
     las=1,
     main="maxPause map3TomaxTSNs KS FDR <= 0.1, e and l >1bp move")
hist(sub_df$AllelicMaxPauseDist,
     freq = F, ylim=c(0,0.4), las=1,
     breaks = seq(-0.5,50,1), col="blue", main="maxPause map3TomaxTSNs KS FDR <= 0.1",
     add=T
)

hist(sub_df_SNP_3G$AllelicMaxPauseDist,     freq = F, ylim=c(0,0.4), las=1,
     breaks = seq(-0.5,50,1),
     col="dark orange", density = 45, add=T)

lines(ecdf(sub_df$AllelicMaxPauseDist), col="blue")
lines(ecdf(sub_df_SNP_3G$AllelicMaxPauseDist), col="dark orange")

legend("right", 
       legend = c( #paste("SNP at early pause, n=", dim(sub_df_SNP)[1], sep=""), 
                   paste("SNP at -3 early pause, n=", dim(sub_df_SNP_3G)[1], sep=""),
                   paste("All FDR<=0.1, >1bp move, n = ", dim(sub_df)[1], sep="")),
       title = ,
       #pch=c(15,15),
       #cex=2, 
       #lty=c(0,0),
       #bty="n",
       lwd=1.5, 
       #density=c(10000,25),
       #angle=c(180,45),
       #angle=45,
       #fill=c("blue","dark orange","dark green")
       col=c("dark orange","blue"),
       pch=16,
       bty = "n"
)


legend("left", 
       legend = c( paste("All FDR<=0.1, >1bp move, n = ", dim(sub_df)[1], sep=""), 
                   paste("SNP at -3 early pause, n=", dim(sub_df_SNP_3G)[1], sep="")),
       title = ,
       #pch=c(15,15),
       #cex=2, 
       lty=c(0,0),
       #bty="n",
       lwd=1.5, 
       density=c(10000,25),
       angle=c(180,45),
       #angle=45,
       fill=c("blue","red")
       , bty = "n"
)
ks.test(sub_df$AllelicMaxPauseDist ,sub_df_SNP_3G$AllelicMaxPauseDist, alternative = "less")



###
# SNPs at -2
# exclude difference 0
snp_position=-2
df$AllelicMaxPauseDist = abs(df$mat_maxPause_map3 - df$pat_maxPause_map3)
sub_df=df[df$map3.p.value.fdr<=0.1 & df$AllelicMaxPauseDist >0,]
sub_df_SNP_2=df[df$map3.p.value.fdr<=0.1 & df$AllelicMaxPauseDist >0 &((df$SNP1_distance_earlyPause==snp_position & !is.na(df$SNP1_distance_earlyPause))|(df$SNP2_distance_earlyPause==snp_position & !is.na(df$SNP2_distance_earlyPause))),]

bed=sub_df_SNP_2

par(mfrow=c(2,1))
plot(ecdf(sub_df$AllelicMaxPauseDist), col="blue", 
     xlab="AllelicMaxPauseDist",
     ylab="density",
     las=1,
     main="maxPause map3TomaxTSNs KS FDR <= 0.1, e and l >1bp move")
hist(sub_df$AllelicMaxPauseDist,
     freq = F, ylim=c(0,0.4), las=1,
     breaks = seq(-0.5,50,1), col="blue", main="maxPause map3TomaxTSNs KS FDR <= 0.1",
     add=T
)

hist(sub_df_SNP_2$AllelicMaxPauseDist,     freq = F, ylim=c(0,0.4), las=1,
     breaks = seq(-0.5,50,1),
     col="dark orange", density = 45, add=T)

lines(ecdf(sub_df$AllelicMaxPauseDist), col="blue")
lines(ecdf(sub_df_SNP_2$AllelicMaxPauseDist), col="dark orange")

legend("right", 
       legend = c( paste("SNP at early pause, n=", dim(sub_df_SNP)[1], sep=""), 
                   paste("SNP at -2 early pause, n=", dim(sub_df_SNP_2)[1], sep=""),
                   paste("All FDR<=0.1, >1bp move, n = ", dim(sub_df)[1], sep="")),
       title = ,
       #pch=c(15,15),
       #cex=2, 
       #lty=c(0,0),
       #bty="n",
       lwd=1.5, 
       #density=c(10000,25),
       #angle=c(180,45),
       #angle=45,
       #fill=c("blue","dark orange","dark green")
       col=c("red","dark orange","blue"),
       pch=16,
       bty = "n"
)


legend("right", 
       legend = c( paste("All FDR<=0.1, >1bp move, n = ", dim(sub_df)[1], sep=""), 
                   paste("SNP at -2 early pause, n=", dim(sub_df_SNP_2)[1], sep="")),
       title = ,
       #pch=c(15,15),
       #cex=2, 
       lty=c(0,0),
       #bty="n",
       lwd=1.5, 
       density=c(10000,25),
       angle=c(180,45),
       #angle=45,
       fill=c("blue","red")
       , bty = "n"
)
ks.test(sub_df$AllelicMaxPauseDist ,sub_df_SNP_2$AllelicMaxPauseDist, alternative = "less")



sum(sub_df$map3.p.value.fdr<=0.1& sub_df$dist_Indel_maxTSN <= sub_df$latePause)
# with indel
sum(sub_df_SNP_3G$map3.p.value.fdr<=0.1& sub_df_SNP_3G$dist_Indel_maxTSN <= sub_df_SNP_3G$latePause)
#View(sub_df_SNP_3G[sub_df_SNP_3G$map3.p.value.fdr<=0.1& sub_df_SNP_3G$dist_Indel_maxTSN <= sub_df_SNP_3G$latePause,])
# without indel
sum(sub_df_SNP_3G$map3.p.value.fdr<=0.1& sub_df_SNP_3G$dist_Indel_maxTSN > sub_df_SNP_3G$latePause)

sum(sub_df_SNP_2$map3.p.value.fdr<=0.1& sub_df_SNP_2$dist_Indel_maxTSN <= sub_df_SNP_2$latePause)
View(sub_df_SNP_2[sum(sub_df_SNP_2$map3.p.value.fdr<=0.1& sub_df_SNP_2$dist_Indel_maxTSN <= sub_df_SNP_2$latePause),])
# without indel
sum(sub_df_SNP_2$map3.p.value.fdr<=0.1& sub_df_SNP_2$dist_Indel_maxTSN > sub_df_SNP_2$latePause)
View(sub_df_SNP_2[(sub_df_SNP_2$map3.p.value.fdr<=0.1& sub_df_SNP_2$dist_Indel_maxTSN > sub_df_SNP_2$latePause),])

# 
SNP_count_and_or <- function(snp_position, subDF=df){
  a = sum(#subDF$map3.p.value.fdr<=0.1&  # KS test 
          #  subDF$AllelicMaxPauseDist >0 & # at least 1bp apart
            ((subDF$SNP1_distance_earlyPause %in% snp_position & !is.na(subDF$SNP1_distance_earlyPause)) 
           # &(subDF$SNP3_distance_earlyPause %in% snp_position & !is.na(subDF$SNP3_distance_earlyPause)) 
              & (subDF$SNP2_distance_earlyPause %in% snp_position & !is.na(subDF$SNP2_distance_earlyPause))))
  o = sum(#subDF$map3.p.value.fdr<=0.1&  # KS test 
          #  subDF$AllelicMaxPauseDist >0 & # at least 1bp apart
            ((subDF$SNP1_distance_earlyPause %in% snp_position & !is.na(subDF$SNP1_distance_earlyPause)) 
             | (subDF$SNP2_distance_earlyPause %in% snp_position & !is.na(subDF$SNP2_distance_earlyPause))))
  count <-NULL
  for (s in snp_position){
    count = c(count, sum(#subDF$map3.p.value.fdr<=0.1&  # KS test 
                        #   subDF$AllelicMaxPauseDist >0 & # at least 1bp apart
                           ((subDF$SNP1_distance_earlyPause==s & !is.na(subDF$SNP1_distance_earlyPause)) 
                             |(subDF$SNP3_distance_earlyPause==s & !is.na(subDF$SNP3_distance_earlyPause)) 
                            | (subDF$SNP2_distance_earlyPause==s & !is.na(subDF$SNP2_distance_earlyPause)))))
  }
  
  return (list(a,o,snp_position, count))
  # return and, or, snp position, count of sites with snp in that position
}


SNP_count_and_or <- function(snp_position, subDF=df){
  a = (subDF[((subDF$SNP1_distance_earlyPause %in% snp_position & !is.na(subDF$SNP1_distance_earlyPause)) 
      #&(subDF$SNP3_distance_earlyPause %in% snp_position & !is.na(subDF$SNP3_distance_earlyPause)) 
     & (subDF$SNP2_distance_earlyPause %in% snp_position & !is.na(subDF$SNP2_distance_earlyPause))),])
  o = 
    (subDF[(subDF$SNP1_distance_earlyPause %in% snp_position & !is.na(subDF$SNP1_distance_earlyPause)) 
                       | (subDF$SNP2_distance_earlyPause %in% snp_position & !is.na(subDF$SNP2_distance_earlyPause))
                       | (subDF$SNP3_distance_earlyPause %in% snp_position & !is.na(subDF$SNP3_distance_earlyPause)),])
  count <-NULL
  for (s in snp_position){
    count = c(count, 
              UniqRowCount(subDF[((subDF$SNP1_distance_earlyPause==s & !is.na(subDF$SNP1_distance_earlyPause)) 
       |(subDF$SNP3_distance_earlyPause==s & !is.na(subDF$SNP3_distance_earlyPause)) 
       | (subDF$SNP2_distance_earlyPause==s & !is.na(subDF$SNP2_distance_earlyPause))),]))
  }
  
  return (list(UniqRowCount(a),UniqRowCount(o),snp_position, count))
  # return and, or, snp position, count of sites with snp in that position
}

df_target=df[df$map3.p.value.fdr<=0.1&  # KS test 
            df$AllelicMaxPauseDist >0,] # at least 1bp apart]

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
SNP_count_and_or(c(0,-2),df_target)
SNP_count_and_or(c(0,-3),df_target)
SNP_count_and_or(c(1:10),df_target)
SNP_count_and_or(c(1,-3),df_target)
SNP_count_and_or(c(-2, 1),df_target)

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
SNP_count_and_or(-3, subDF = df_target[df_target$dist_Indel_maxTSN > df_target$latePause,])
SNP_count_and_or(1, subDF = df_target[df_target$dist_Indel_maxTSN > df_target$latePause,])
# without indel AND with any snps
subDF=df_target[df_target$dist_Indel_maxTSN > df_target$latePause,]
sum(!is.na(subDF$SNP1_distance_earlyPause))

# pause distribution not significantly different
df_target=df[df$map3.p.value.fdr> 0.9 & # KS test 
               df$AllelicMaxPauseDist ==0,] # at least 1bp apart]
sum(df$map3.p.value.fdr> 0.9)
dim(df_target)

# SNPs
SNP_count_and_or(c(0,-2),df_target) 
SNP_count_and_or(c(0,-3),df_target)
SNP_count_and_or(c(-2, -3),df_target)
#SNP_count_and_or(c(0,-2, -3),df_target)


# with indel
sum(df_target$dist_Indel_maxTSN <= df_target$latePause)
SNP_count_and_or(0, subDF = df_target[df_target$dist_Indel_maxTSN <= df_target$latePause,])
SNP_count_and_or(-2, subDF = df_target[df_target$dist_Indel_maxTSN <= df_target$latePause,])
SNP_count_and_or(-3, subDF = df_target[df_target$dist_Indel_maxTSN <= df_target$latePause,])
#SNP_count_and_or(c(0,-2, -3),subDF = df_target[df_target$dist_Indel_maxTSN <= df_target$latePause,])
# with indel and any snps
subDF=df_target[df_target$dist_Indel_maxTSN <= df_target$latePause,]
sum(!is.na(subDF$SNP1_distance_earlyPause))


# without indel
SNP_count_and_or(0, subDF = df_target[df_target$dist_Indel_maxTSN > df_target$latePause,])
SNP_count_and_or(-2, subDF = df_target[df_target$dist_Indel_maxTSN > df_target$latePause,])
SNP_count_and_or(-3, subDF = df_target[df_target$dist_Indel_maxTSN > df_target$latePause,])
#SNP_count_and_or(c(0,-2, -3),subDF = df_target[df_target$dist_Indel_maxTSN > df_target$latePause,])

# without indel AND with any snps other than 0,-2,-3
subDF=df_target[df_target$dist_Indel_maxTSN > df_target$latePause,]
sum(!is.na(subDF$SNP1_distance_earlyPause))


# significant sites without SNPs at 0,-2,-3
df_target=df[df$map3.p.value.fdr<=0.1&  # KS test 
               df$AllelicMaxPauseDist >0,] # at least 1bp apart]
# no indel
df_target = df_target[df_target$dist_Indel_maxTSN > df_target$latePause,]
# no snps at 0,-2,-3
df_target =df_target[!df_target$SNP1_distance_earlyPause %in% c(0,-2,-3),] 
dim(df_target)
df_target =df_target[!df_target$SNP2_distance_earlyPause %in% c(0,-2,-3),] 
dim(df_target)
df_target =df_target[!df_target$SNP3_distance_earlyPause %in% c(0,-2,-3),] 
dim(df_target)
View(df_target)  

## SNPs around early pause of significant sites without indel and SNPs at 0,-2,-3
par(mfrow=c(3,1))
step=1
metaplot.SNPsLocation.aroundEarlyPause(df_target, show.window=100,
                                       name=paste("HSK, step=", step,", early and late at least ", bp.apart," bp apart, upto " ,upto, "no indel, no 0,-2,-3 SNPs", sep=""), col="red", step=step)
abline(v=0, col="gray")

metaplot.SNPsLocation.aroundEarlyPause(df[df$map3.p.value.fdr>0.9,], 
                                       show.window=100,
                                       col="blue", step=step, add=TRUE)

legend("topright", legend=c(paste("Early Pause, subset of fdr <= 0.1, n= ", dim(df_target)[1], sep=""), 
                            paste("Early Pause, fdr >  0.9, n= ", sum(df$map3.p.value.fdr>0.9), sep="")),
       col=c("red", "blue"), 
       bty = "n", lty=1, pch=19)

# add a test
e1=metaplot.SNPsLocation.aroundEarlyPause(df_target, show.window=100,
                                       name=paste("HSK, step=", step,", early and late at least ", bp.apart," bp apart, upto " ,upto, "no indel, no 0,-2,-3 SNPs", sep=""), col="red", step=step,
                                       use.sum = TRUE, plot = FALSE)
e9=metaplot.SNPsLocation.aroundEarlyPause(df[df$map3.p.value.fdr>0.9,], col="blue", step=step, 
                                          show.window=100,
                                          use.sum = TRUE, plot = FALSE)

fdr_cutoff=0.1
p.value <- NULL
odds.ratio <- NULL
testors <- NULL

for (n in (1:length(e1[[2]]))){
  
  testor <-  matrix(c(e1[[2]][n], dim(df_target)[1],
                      e9[[2]][n], sum(df$map3.p.value.fdr>0.9)),
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
# text(e1[[1]][which(p.value <= 0.1)],p.value[which(p.value <= 0.1)], label=paste(round(p.value[which(p.value <= 0.1)], digits = 2), sep=" "))
# dev.off()