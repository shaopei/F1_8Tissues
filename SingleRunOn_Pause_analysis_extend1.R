
t="HT"
show.window=100
pause_window_0.1 <- read.table(paste("/Volumes/SPC_SD/KD_IGV/",t,"_dREG_5mat5pat_uniq_pValue_fdr0.1.bed.bed", sep =""), header = F)
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
  
  t=c(map5.pat.peak[i],   map5.mat.peak[i])   # the base with of max map5 reads 
  
  # determine which allele is early pause, which is late pause
  a=sort.int(c(hmat.pat.peak[i],  hmat.mat.peak[i]), index.return = TRUE )
  AT$early.pause[i] = a$x[1] 
  AT$late.pause[i] = a$x[2]
  AT$TSN.early.pause[i] = t[a$ix[1]]  # from the same allele as early pause, which base has the max map5 reads
  AT$TSN.late.pause[i] = t[a$ix[2]]   # from the same allele as late pause, which base has the max map5 reads,
}
### from heatmap.Pause function end ###


# distribution of distance between early and late TSN
AT$early.TSN.pause.dist = AT$early.pause - AT$TSN.early.pause
AT$late.TSN.pause.dist = AT$late.pause - AT$TSN.late.pause
AT$valid.pairs = AT$early.TSN.pause.dist > 0 & AT$late.TSN.pause.dist > 0
old_AT = AT
AT = AT[AT$valid.pairs,]

AT$TSN.dist = AT$TSN.late.pause - AT$TSN.early.pause
View(AT)
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
hist(AT$pause.dist[AT$pause.dist!= 0], main="distance between early and late pause, exclude dist==0")
sum(AT$pause.dist==0)
sum(AT$pause.dist!=0)
sum(AT$pause.dist> 10)

plot(AT$TSN.dist, AT$pause.dist)
plot(abs(AT$TSN.dist), AT$pause.dist, xlim=c(0,100), ylim=c(0,100), pch=19, col=rgb(0,0,0,alpha = 0.25))
abline(0,1, col="blue")
plot(AT$TSN.dist, AT$pause.dist, xlim=c(0,50), ylim=c(0,50), pch=19, col=rgb(0,0,0,alpha = 0.125))
abline(0,1)

# modified "AT$valid.pairs" in heatmap.Pause function
# AT$valid.pairs = (AT$early.TSN.pause.dist > 0) & (AT$late.TSN.pause.dist > 0)

t="CenteratTSN_sortBydistEarlyLatePause_HT_fdr0.1_TSN.group=T"
heatmap.Pause_allreads(pause_window_0.1, file.plus.bw, file.minus.bw ,
                       dist=200, step=2, up_dist =100,
                       file.pdf=paste(t,"_allreads_heatmap_step2_up100_dist200_map5.pdf",sep = ""),
                       metaplot=FALSE,
                       metaplot.pdf = paste(t,"_allreads_metaplot_step2_up100_dist200_map5.pdf",sep = ""),
                       allelic.heatmap.pdf = paste(t,"_allelicReads_heatmap_step2_up100_dist200_map5.pdf",sep = ""),
                       allelic.metaplot.pdf = paste(t,"_allelicReads__metaplot_step2_up100_dist200_map5.pdf",sep = ""),
                       breaks=seq(0, 20, 0.001), times=times, map5=TRUE, metaplot.ylim = c(0,160),
                       TSN.group=TRUE, center.at.TSN = TRUE)

t="CenteratTSN_sortBydistEarlyLatePause_HT_fdr0.1_TSN.group=FALSE"
heatmap.Pause_allreads(pause_window_0.1, file.plus.bw, file.minus.bw ,
                       dist=200, step=2, up_dist =100,
                       file.pdf=paste(t,"_allreads_heatmap_step2_up100_dist200_map5.pdf",sep = ""),
                       metaplot.pdf = paste(t,"_allreads_metaplot_step2_up100_dist200_map5.pdf",sep = ""),
                       metaplot=FALSE,
                       allelic.heatmap.pdf = paste(t,"_allelicReads_heatmap_step2_up100_dist200_map5.pdf",sep = ""),
                       allelic.metaplot.pdf = paste(t,"_allelicReads__metaplot_step2_up100_dist200_map5.pdf",sep = ""),
                       breaks=seq(0, 20, 0.001), times=times, map5=TRUE, metaplot.ylim = c(0,160),
                       TSN.group=FALSE, center.at.TSN = TRUE)


t="CenteratEarlyPause_sortBydistEarlyLatePause_HT_fdr0.1_TSN.group=T"
heatmap.Pause_allreads(pause_window_0.1, file.plus.bw, file.minus.bw ,
                       dist=200, step=2, up_dist =100,
                       file.pdf=paste(t,"_allreads_heatmap_step2_up100_dist200_map5.pdf",sep = ""),
                       metaplot=FALSE,
                       metaplot.pdf = paste(t,"_allreads_metaplot_step2_up100_dist200_map5.pdf",sep = ""),
                       allelic.heatmap.pdf = paste(t,"_allelicReads_heatmap_step2_up100_dist200_map5.pdf",sep = ""),
                       allelic.metaplot.pdf = paste(t,"_allelicReads__metaplot_step2_up100_dist200_map5.pdf",sep = ""),
                       breaks=seq(0, 20, 0.001), times=times, map5=TRUE, metaplot.ylim = c(0,160),
                       TSN.group=TRUE, center.at.TSN = FALSE)

t="CenteratEarlyPause_sortBydistEarlyLatePause_HT_fdr0.1_TSN.group=FALSE"
heatmap.Pause_allreads(pause_window_0.1, file.plus.bw, file.minus.bw ,
                       dist=200, step=2, up_dist =100,
                       file.pdf=paste(t,"_allreads_heatmap_step2_up100_dist200_map5.pdf",sep = ""),
                       metaplot=FALSE,
                       metaplot.pdf = paste(t,"_allreads_metaplot_step2_up100_dist200_map5.pdf",sep = ""),
                       allelic.heatmap.pdf = paste(t,"_allelicReads_heatmap_step2_up100_dist200_map5.pdf",sep = ""),
                       allelic.metaplot.pdf = paste(t,"_allelicReads__metaplot_step2_up100_dist200_map5.pdf",sep = ""),
                       breaks=seq(0, 20, 0.001), times=times, map5=TRUE, metaplot.ylim = c(0,160),
                       TSN.group=FALSE, center.at.TSN = FALSE)
###
# modified "AT$valid.pairs" in heatmap.Pause function
# AT$valid.pairs = (AT$early.TSN.pause.dist > 0) & (AT$late.TSN.pause.dist > 0) & (AT$TSN.late.pause - AT$TSN.early.pause > 0)

t="CenteratTSN_sortBydistEarlyLatePause_HT_fdr0.1_TSN.group=T_dist.TSN>0"
heatmap.Pause_allreads(pause_window_0.1, file.plus.bw, file.minus.bw ,
                       dist=200, step=2, up_dist =100,
                       file.pdf=paste(t,"_allreads_heatmap_step2_up100_dist200_map5.pdf",sep = ""),
                       metaplot=FALSE,
                       metaplot.pdf = paste(t,"_allreads_metaplot_step2_up100_dist200_map5.pdf",sep = ""),
                       allelic.heatmap.pdf = paste(t,"_allelicReads_heatmap_step2_up100_dist200_map5.pdf",sep = ""),
                       allelic.metaplot.pdf = paste(t,"_allelicReads__metaplot_step2_up100_dist200_map5.pdf",sep = ""),
                       breaks=seq(0, 20, 0.001), times=times, map5=TRUE, metaplot.ylim = c(0,160),
                       TSN.group=TRUE, center.at.TSN = TRUE)

t="CenteratEarlyPause_sortBydistEarlyLatePause_HT_fdr0.1_TSN.group=T_dist.TSN>0"
heatmap.Pause_allreads(pause_window_0.1, file.plus.bw, file.minus.bw ,
                       dist=200, step=2, up_dist =100,
                       file.pdf=paste(t,"_allreads_heatmap_step2_up100_dist200_map5.pdf",sep = ""),
                       metaplot=FALSE,
                       metaplot.pdf = paste(t,"_allreads_metaplot_step2_up100_dist200_map5.pdf",sep = ""),
                       allelic.heatmap.pdf = paste(t,"_allelicReads_heatmap_step2_up100_dist200_map5.pdf",sep = ""),
                       allelic.metaplot.pdf = paste(t,"_allelicReads__metaplot_step2_up100_dist200_map5.pdf",sep = ""),
                       breaks=seq(0, 20, 0.001), times=times, map5=TRUE, metaplot.ylim = c(0,160),
                       TSN.group=TRUE, center.at.TSN = FALSE)

# modified "AT$valid.pairs" in heatmap.Pause function
# AT$valid.pairs = (AT$early.TSN.pause.dist > 0) & (AT$late.TSN.pause.dist > 0) & (AT$TSN.late.pause == AT$TSN.early.pause)

t="CenteratTSN_sortBydistEarlyLatePause_HT_fdr0.1_TSN.group=T_dist.TSN=0"
heatmap.Pause_allreads(pause_window_0.1, file.plus.bw, file.minus.bw ,
                       dist=200, step=2, up_dist =100,
                       file.pdf=paste(t,"_allreads_heatmap_step2_up100_dist200_map5.pdf",sep = ""),
                       metaplot=FALSE,
                       metaplot.pdf = paste(t,"_allreads_metaplot_step2_up100_dist200_map5.pdf",sep = ""),
                       allelic.heatmap.pdf = paste(t,"_allelicReads_heatmap_step2_up100_dist200_map5.pdf",sep = ""),
                       allelic.metaplot.pdf = paste(t,"_allelicReads__metaplot_step2_up100_dist200_map5.pdf",sep = ""),
                       breaks=seq(0, 20, 0.001), times=times, map5=TRUE, metaplot.ylim = c(0,160),
                       TSN.group=TRUE, center.at.TSN = TRUE)

t="CenteratEarlyPause_sortBydistEarlyLatePause_HT_fdr0.1_TSN.group=T_dist.TSN=0"
heatmap.Pause_allreads(pause_window_0.1, file.plus.bw, file.minus.bw ,
                       dist=200, step=2, up_dist =100,
                       file.pdf=paste(t,"_allreads_heatmap_step2_up100_dist200_map5.pdf",sep = ""),
                       metaplot=FALSE,
                       metaplot.pdf = paste(t,"_allreads_metaplot_step2_up100_dist200_map5.pdf",sep = ""),
                       allelic.heatmap.pdf = paste(t,"_allelicReads_heatmap_step2_up100_dist200_map5.pdf",sep = ""),
                       allelic.metaplot.pdf = paste(t,"_allelicReads__metaplot_step2_up100_dist200_map5.pdf",sep = ""),
                       breaks=seq(0, 20, 0.001), times=times, map5=TRUE, metaplot.ylim = c(0,160),
                       TSN.group=TRUE, center.at.TSN = FALSE)

# modified "AT$valid.pairs" in heatmap.Pause function
#AT$valid.pairs = (AT$early.TSN.pause.dist > 0) & (AT$late.TSN.pause.dist > 0) & (AT$TSN.late.pause == AT$TSN.early.pause) & (AT$early.pause == AT$late.pause)
t="CenteratTSN_sortBydistEarlyLatePause_HT_fdr0.1_TSN.group=T_dist.TSN=0_dist.Pause=0"
heatmap.Pause_allreads(pause_window_0.1, file.plus.bw, file.minus.bw ,
                       dist=200, step=2, up_dist =100,
                       file.pdf=paste(t,"_allreads_heatmap_step2_up100_dist200_map5.pdf",sep = ""),
                       metaplot=FALSE,
                       metaplot.pdf = paste(t,"_allreads_metaplot_step2_up100_dist200_map5.pdf",sep = ""),
                       allelic.heatmap.pdf = paste(t,"_allelicReads_heatmap_step2_up100_dist200_map5.pdf",sep = ""),
                       allelic.metaplot.pdf = paste(t,"_allelicReads__metaplot_step2_up100_dist200_map5.pdf",sep = ""),
                       breaks=seq(0, 20, 0.001), times=times, map5=TRUE, metaplot.ylim = c(0,160),
                       TSN.group=TRUE, center.at.TSN = TRUE)

t="CenteratEarlyPause_sortBydistEarlyLatePause_HT_fdr0.1_TSN.group=T_dist.TSN=0_dist.Pause=0"
heatmap.Pause_allreads(pause_window_0.1, file.plus.bw, file.minus.bw ,
                       dist=200, step=2, up_dist =100,
                       file.pdf=paste(t,"_allreads_heatmap_step2_up100_dist200_map5.pdf",sep = ""),
                       metaplot=FALSE,
                       metaplot.pdf = paste(t,"_allreads_metaplot_step2_up100_dist200_map5.pdf",sep = ""),
                       allelic.heatmap.pdf = paste(t,"_allelicReads_heatmap_step2_up100_dist200_map5.pdf",sep = ""),
                       allelic.metaplot.pdf = paste(t,"_allelicReads__metaplot_step2_up100_dist200_map5.pdf",sep = ""),
                       breaks=seq(0, 20, 0.001), times=times, map5=TRUE, metaplot.ylim = c(0,160),
                       TSN.group=TRUE, center.at.TSN = FALSE)
#    AT$valid.pairs = (AT$early.TSN.pause.dist > 0) & (AT$late.TSN.pause.dist > 0) & (AT$TSN.late.pause == AT$TSN.early.pause) & (AT$early.pause != AT$late.pause)
t="CenteratTSN_sortBydistEarlyLatePause_HT_fdr0.1_TSN.group=T_dist.TSN=0_dist.Pause!=0"
heatmap.Pause_allreads(pause_window_0.1, file.plus.bw, file.minus.bw ,
                       dist=200, step=2, up_dist =100,
                       file.pdf=paste(t,"_allreads_heatmap_step2_up100_dist200_map5.pdf",sep = ""),
                       metaplot=FALSE,
                       metaplot.pdf = paste(t,"_allreads_metaplot_step2_up100_dist200_map5.pdf",sep = ""),
                       allelic.heatmap.pdf = paste(t,"_allelicReads_heatmap_step2_up100_dist200_map5.pdf",sep = ""),
                       allelic.metaplot.pdf = paste(t,"_allelicReads__metaplot_step2_up100_dist200_map5.pdf",sep = ""),
                       breaks=seq(0, 20, 0.001), times=times, map5=TRUE, metaplot.ylim = c(0,160),
                       TSN.group=TRUE, center.at.TSN = TRUE)

t="CenteratEarlyPause_sortBydistEarlyLatePause_HT_fdr0.1_TSN.group=T_dist.TSN=0_dist.Pause!=0"
heatmap.Pause_allreads(pause_window_0.1, file.plus.bw, file.minus.bw ,
                       dist=200, step=2, up_dist =100,
                       file.pdf=paste(t,"_allreads_heatmap_step2_up100_dist200_map5.pdf",sep = ""),
                       metaplot=FALSE,
                       metaplot.pdf = paste(t,"_allreads_metaplot_step2_up100_dist200_map5.pdf",sep = ""),
                       allelic.heatmap.pdf = paste(t,"_allelicReads_heatmap_step2_up100_dist200_map5.pdf",sep = ""),
                       allelic.metaplot.pdf = paste(t,"_allelicReads__metaplot_step2_up100_dist200_map5.pdf",sep = ""),
                       breaks=seq(0, 20, 0.001), times=times, map5=TRUE, metaplot.ylim = c(0,160),
                       TSN.group=TRUE, center.at.TSN = FALSE)
