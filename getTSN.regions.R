# try to identify the regions showing allelic difference in initiation distibution

setwd("~/Box Sync/Danko_lab_work/F1_8Tissues/Kidney_and_SingleRunOn/Pause_manuscript")
source("/Users/shaopei/Box\ Sync/Danko_lab_work/F1_8Tissues/PolyA_Allele-specific/heatmap/heatmaps.R")


t="HT"
pause_window_0.1 <- read.table(paste("/Volumes/SPC_SD/KD_IGV/TSS/",t,"_dREG_5mat5pat_uniq_pValue_fdr0.1_TSS.bed", sep =""), header = F)
pause_window_0.1 <- pause_window_0.1[pause_window_0.1$V1 != 'chrX',]
end=".rpm.bw"; times=10
#end=".bw"; times=1
# allelic reads (map3)
file.bw.plus.pat <- paste("/Volumes/SPC_SD/KD_IGV/",t,".pat.map2ref.1bp_plus",end, sep="")
file.bw.minus.pat <- paste("/Volumes/SPC_SD/KD_IGV/",t,".pat.map2ref.1bp_minus",end, sep="")
file.bw.plus.mat <- paste("/Volumes/SPC_SD/KD_IGV/",t,".mat.map2ref.1bp_plus",end, sep="")
file.bw.minus.mat <- paste("/Volumes/SPC_SD/KD_IGV/",t,".mat.map2ref.1bp_minus",end, sep="")
# allelic reads (map3) HT.mat.map5.map2ref.1bp_minus.bw
map5.file.bw.plus.pat <- paste("/Volumes/SPC_SD/KD_IGV/map5/",t,".pat.map5.map2ref.1bp_plus.bw", sep="")
map5.file.bw.minus.pat <- paste("/Volumes/SPC_SD/KD_IGV/map5/",t,".pat.map5.map2ref.1bp_minus.bw", sep="")
map5.file.bw.plus.mat <- paste("/Volumes/SPC_SD/KD_IGV/map5/",t,".mat.map5.map2ref.1bp_plus.bw", sep="")
map5.file.bw.minus.mat <- paste("/Volumes/SPC_SD/KD_IGV/map5/",t,".mat.map5.map2ref.1bp_minus.bw", sep="")
# all reads
file.plus.bw <- paste("/Volumes/SPC_SD/KD_IGV/",t,"_PB6_F5N6_dedup_QC_end_plus",end, sep="")
file.minus.bw <- paste("/Volumes/SPC_SD/KD_IGV/",t,"_PB6_F5N6_dedup_QC_end_minus",end, sep="")
map5.file.plus.bw <- paste("/Volumes/SPC_SD/KD_IGV/map5/",t,"_PB6_F5N6_dedup_QC_end_map5_plus.bw", sep="")
map5.file.minus.bw <- paste("/Volumes/SPC_SD/KD_IGV/map5/",t,"_PB6_F5N6_dedup_QC_end_map5_minus.bw", sep="")
SNP.bw <- "/Volumes/SPC_SD/KD_IGV/P.CAST_M.B6_indelsNsnps_CAST.bam.snp.unfiltered_plus.bw"
AT=pause_window_0.1


AT <- AT[,1:6] 

# identify the location of max reads in mat and pat pause
# determine the location of early and late TSN

map5.pat.peak <- NULL
map5.mat.peak <- NULL
hmat.pat.peak <- NULL  #map3
hmat.mat.peak <- NULL  # map3
for (i in 1:NROW(AT)){
  # identify the base with of max map5 reads (maxTSN)
  map5.pat <- read_read_mat_S (map5.file.bw.plus.pat, map5.file.bw.minus.pat, AT[i,], step, times=times, use.log=use.log)
  map5.mat <- read_read_mat_S (map5.file.bw.plus.mat, map5.file.bw.minus.mat, AT[i,], step, times=times, use.log=use.log) 
  map5.pat.peak[i] <- which(map5.pat == max(map5.pat))  # if two bases with same read count, output the 5 prime one (up stream)
  map5.mat.peak[i] <- which(map5.mat == max(map5.mat)) 
  # identify the base with of max map3 reads  (maxPause)
  hmat.pat <- read_read_mat_S (file.bw.plus.pat, file.bw.minus.pat, AT[i,], step, times=times, use.log=use.log)
  hmat.mat <- read_read_mat_S (file.bw.plus.mat, file.bw.minus.mat, AT[i,], step, times=times, use.log=use.log) 
  hmat.pat.peak[i] <- which(hmat.pat == max(hmat.pat))
  hmat.mat.peak[i] <- which(hmat.mat == max(hmat.mat))
  t=c(hmat.pat.peak[i],  hmat.mat.peak[i])   # the base with of max map3 reads 
  # determine which allele is early TSN, which is late TSN
  a=sort.int(c(map5.pat.peak[i],   map5.mat.peak[i]), index.return = TRUE )
  AT$early.TSN[i] = a$x[1] 
  AT$late.TSN[i] = a$x[2]
  AT$early.TSN.pause[i] = t[a$ix[1]]  # from the same allele as early TSN, which base has the max map3 reads
  AT$late.TSN.pause[i] = t[a$ix[2]]   # from the same allele as late TSN, which base has the max map3 reads,
}


# distribution of distance between early and late TSN
AT$early.TSN.pause.dist = AT$early.TSN.pause - AT$early.TSN
AT$late.TSN.pause.dist = AT$late.TSN.pause - AT$late.TSN
AT$valid.pairs = AT$early.TSN.pause.dist > 0 & AT$late.TSN.pause.dist > 0
dim(AT)
sum(AT$valid.pairs)

AT$TSN.dist = AT$late.TSN - AT$early.TSN
View(AT)
hist(AT$TSN.dist, main="distance between early and late TSN")
hist(AT$TSN.dist[AT$TSN.dist !=0 ], main="distance between early and late TSN, exclude dist==0")
sum(AT$TSN.dist == 0)
sum(AT$TSN.dist !=0 )
sum(AT$TSN.dist <0 )
sum(AT$TSN.dist >0 )
sum(AT$TSN.dist >10 )


old_AT = AT
AT = AT[AT$valid.pairs,]

# distribution of distance between early and late pause
AT$pause.dist= AT$late.TSN.pause - AT$early.TSN.pause
hist(AT$pause.dist, main="distance between early and late pause")
hist(AT$pause.dist[AT$pause.dist!= 0], main="distance between early and late pause, exclude dist==0")
sum(AT$pause.dist==0)
sum(AT$pause.dist < 0)
sum(AT$pause.dist> 0)
sum(AT$pause.dist> 10)

plot(AT$TSN.dist, AT$pause.dist)
plot(AT$TSN.dist, (AT$pause.dist), xlim=c(0,100), ylim=c(-100,100), pch=19, col=rgb(0,0,0,alpha = 0.25))
abline(0,1, col="blue")
plot(AT$TSN.dist, AT$pause.dist, xlim=c(0,50), ylim=c(0,50), pch=19, col=rgb(0,0,0,alpha = 0.125))
abline(0,1)

