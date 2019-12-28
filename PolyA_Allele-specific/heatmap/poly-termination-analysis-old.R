setwd("~/Box Sync/Danko_lab_work/F1_8Tissues/PolyA_Allele-specific/heatmap/")
library(bigWig)

read_count <- function(file.plus.bw, file.minus.bw , bed6) {
    bw.plus  <- load.bigWig( file.plus.bw )
    bw.minus <- load.bigWig( file.minus.bw )
    CountMatrix <- bed6.region.bpQuery.bigWig(bw.plus, bw.minus, bed6[,c(1:6)] ,abs.value=TRUE, op = "sum")
    unload.bigWig(bw.plus);
    unload.bigWig(bw.minus);
    
    return(CountMatrix);    
}

Tunit_BW_pairs <- function (at, t){
tunit <-  read.table(paste("../",at,"_AT_4tunitIntersectNativeHMM_tunits.bed", sep=""), header = F)
file.plus.bw <- paste("../",t,"_all_plus.rpm.bw", sep="")
file.minus.bw <- paste("../",t,"_all_minus.rpm.bw", sep="")
log10(read_count(file.plus.bw, file.minus.bw, tunit)*1000/(tunit$V3 - tunit$V2)+1) # log(rpkm+1)
}

BN.Tunit_BN.bw <- Tunit_BW_pairs ("BN","BN")
BN.Tunit_LV.bw <- Tunit_BW_pairs ("BN","LV")
LV.Tunit_LV.bw <- Tunit_BW_pairs ("LV","LV")
LV.Tunit_BN.bw <- Tunit_BW_pairs ("LV","BN")
summary(BN.Tunit_BN.bw)

library("vioplot")
library(RColorBrewer)
#display.brewer.all()
#display.brewer.pal(n = 9,name = "PuOr")
PuOr = brewer.pal(n = 9,name = "PuOr")


BN_col = PuOr[7]
LV_col = PuOr[3] #'#E69F00'
BN_l = PuOr[8]
LV_l = PuOr[2] #'#E69F00'
vioplot(BN.Tunit_BN.bw, BN.Tunit_LV.bw, LV.Tunit_BN.bw, LV.Tunit_LV.bw,
        names=c("BN.BN", "BN.LV", "LV.BN", "LV.LV"), 
        border=c(BN_l, BN_l, LV_l, LV_l), lwd=2,
        col=c(BN_col, LV_col,BN_col, LV_col),
        xlab="Tunit.BW",
        ylab="log10(rpkm+1)",
        #ylim = c(-1,9),
        las=2)

legend("topleft", legend=c("Brain", "Liver"),
       fill=c(BN_col, LV_col), cex=1.5 , bty = "n",
       border = c(BN_l, LV_l))


# vioplot( LV.Tunit_BN.bw, LV.Tunit_LV.bw,BN.Tunit_BN.bw, BN.Tunit_LV.bw,
#         names=c( "LV.BN", "LV.LV", "BN.BN", "BN.LV"),
#         border=c(LV_l, LV_l,BN_l, BN_l), lwd=2,
#         col=c(BN_col, LV_col,BN_col, LV_col),
#         xlab="Tunit.BW",
#         ylab="log(rpkm+1)",
#         #ylim = c(-1,9),
#         las=2)


Tunit_BW_pairs_rpkm <- function (at, t){
  tunit <-  read.table(paste("../",at,"_AT_4tunitIntersectNativeHMM_tunits.bed", sep=""), header = F)
  file.plus.bw <- paste("../",t,"_all_plus.rpm.bw", sep="")
  file.minus.bw <- paste("../",t,"_all_minus.rpm.bw", sep="")
  read_count(file.plus.bw, file.minus.bw, tunit)*1000/(tunit$V3 - tunit$V2) 
}

BN.Tunit_BN.bw <- Tunit_BW_pairs_rpkm  ("BN","BN")
BN.Tunit_LV.bw <- Tunit_BW_pairs_rpkm  ("BN","LV")
LV.Tunit_LV.bw <- Tunit_BW_pairs_rpkm  ("LV","LV")
LV.Tunit_BN.bw <- Tunit_BW_pairs_rpkm  ("LV","BN")

par(mfrow=c(2,1))
hist(log2(BN.Tunit_LV.bw / BN.Tunit_BN.bw), xlim=c(-12,10), freq = F, ylim=c(0,0.4))
hist(log2(LV.Tunit_BN.bw / LV.Tunit_LV.bw), xlim=c(-12,10), freq = F, ylim=c(0,0.4))

par(mfrow=c(4,1))
hist(log10(BN.Tunit_BN.bw+0.01), breaks = seq(-2,2,0.1), border = BN_l, col = BN_col, ylim = c(0,100))
hist(log10(BN.Tunit_LV.bw+0.01), breaks = seq(-2,2,0.1), border = BN_l, col=LV_col, ylim = c(0,100))
hist(log10(LV.Tunit_BN.bw+0.01), breaks = seq(-2,2,0.1), border = LV_l, col = BN_col, ylim = c(0,100))
hist(log10(LV.Tunit_LV.bw+0.01), breaks = seq(-2,2,0.1), border = LV_l, col=LV_col, ylim = c(0,100))

dev.off()
vioplot(log10(BN.Tunit_BN.bw+0.001),
        log10(BN.Tunit_LV.bw+0.001),
        log10(LV.Tunit_BN.bw+0.001),
        log10(LV.Tunit_LV.bw+0.001),
        las=1,
        #horizontal =T, 
        names = c("log10(BN.Tunit_BN.bw)",
                                 "log10(BN.Tunit_LV.bw)",
                                 "log10(LV.Tunit_BN.bw)",
                                 "log10(LV.Tunit_LV.bw+0.001)"),
        border=c(BN_l, BN_l, LV_l, LV_l), lwd=2,
        col=c(BN_col, LV_col,BN_col, LV_col), 
        ylim = c(-3,2))

## annotate each tunit table
navg=1
#at="LV"; other = "BN"
at="BN"; other = "LV"
AT <-  read.table(paste("../",at,"_AT_3tunitIntersectNativeHMM_intersectRegion.bed", sep=""), header = F)
AT$value <-Tunit_BW_pairs_rpkm  (at,other) /  Tunit_BW_pairs_rpkm  (at,at)
AT <- AT[order(AT$value, decreasing = TRUE),]

for (t in c("LV", "BN")){
  file.bw.plus.pat <- paste("../",t,"_MB6_all_R1.pat_1bp_plus.bw", sep="")
  file.bw.minus.pat <- paste("../",t,"_MB6_all_R1.pat_1bp_minus.bw", sep="")
  file.bw.plus.mat <- paste("../",t,"_MB6_all_R1.mat_1bp_plus.bw", sep="")
  file.bw.minus.mat <- paste("../",t,"_MB6_all_R1.mat_1bp_minus.bw", sep="")
  hl.bw.plus.pat <- paste("../",at,"_MB6_all_R1.pat_1bp_plus.bw", sep="")
  hl.bw.minus.pat <- paste("../",at,"_MB6_all_R1.pat_1bp_minus.bw", sep="")
  hl.bw.plus.mat <- paste("../",at,"_MB6_all_R1.mat_1bp_plus.bw", sep="")
  hl.bw.minus.mat <- paste("../",at,"_MB6_all_R1.mat_1bp_minus.bw", sep="")
  #file.plus.bw <- paste("../",t,"_all_plus.rpm.bw", sep="")
  #file.minus.bw <- paste("../",t,"_all_minus.rpm.bw", sep="")
 
  s = seq(1,dim(AT)[1], dim(AT)[1]%/%4)

for (i in 1:4 ){
  q <- AT[s[i]:s[i+1],]
  #heatmap.AT3(q, file.bw.plus.pat,file.bw.minus.pat, file.bw.plus.mat ,file.bw.minus.mat, dist=50000, step=500,  file.pdf=paste("BN","-ATtunitAlleleHMMInterectRegion_", t,"-bw-heatmap-allReads_AT50Kb_step500_navg","10","_q",i,".pdf",sep=""), bl_wd=1, show.AT.line=TRUE, navg=navg, use.log=FALSE, times=10, breaks = seq(0,10,0.1))
  heatmap.AT3(q, file.bw.plus.pat,file.bw.minus.pat, file.bw.plus.mat ,file.bw.minus.mat, hl.bw.plus.pat=hl.bw.plus.pat ,hl.bw.minus.pat=hl.bw.minus.pat, hl.bw.plus.mat = hl.bw.plus.mat ,hl.bw.minus.mat=hl.bw.minus.mat, dist=50000, step=500,  file.pdf=paste(t,"_q",i,".pdf",sep=""), bl_wd=1, show.AT.line=TRUE, navg=navg, use.log=FALSE, times=1, breaks = seq(0,10/i,0.01))
  }
}


### the relationship between Tunit expression level and AT window length
setwd("~/Box Sync/Danko_lab_work/F1_8Tissues/PolyA_Allele-specific/heatmap/")
at = "LV"
tunit <-  read.table(paste("../",at,"_AT_3tunitIntersectNativeHMM_tunits.bed", sep=""), header = F)
AT <- read.table(paste("../",at,"_AT_3tunitIntersectNativeHMM_intersectRegion.bed", sep=""), header = F)
file.plus.bw <- paste("../",at,"_all_plus.rpm.bw", sep="")
file.minus.bw <- paste("../",at,"_all_minus.rpm.bw", sep="")
tunit_rpkm <- read_count(file.plus.bw, file.minus.bw, tunit) *1000/(tunit$V3 - tunit$V2) 
tunit_rpmSum <- read_count(file.plus.bw, file.minus.bw, tunit) 
 smoothScatter(log10(AT$V3-AT$V2), log(tunit_rpkm))
 points(log10(AT$V3-AT$V2), log(tunit_rpkm))
# 
# smoothScatter(log(AT$V3-AT$V2), log(tunit_rpmSum))
# points(log(AT$V3-AT$V2), log(tunit_rpmSum))
# 
# plot_colorByDensity = function(x1,x2,
#                                ylim=c(min(x2),max(x2)),
#                                xlim=c(min(x1),max(x1)),
#                                xlab="",ylab="",main="") {
#   
#   df <- data.frame(x1,x2)
#   x <- densCols(x1,x2, colramp=colorRampPalette(c("black", "white")))
#   df$dens <- col2rgb(x)[1,] + 1L
#   cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F","#FCFF00", "#FF9400", "#FF3100"))(256)
#   df$col <- cols[df$dens]
#   plot(x2~x1, data=df[order(df$dens),], 
#        ylim=ylim,xlim=xlim,pch=20,col=col,
#        cex=2,xlab=xlab,ylab=ylab,
#        main=main)
# }
# 
# plot_colorByDensity(log10(AT$V3-AT$V2), log(tunit_rpkm))
# smoothScatter(log10(AT$V3-AT$V2), log(tunit_rpkm))
# points(log10(AT$V3-AT$V2), log(tunit_rpkm))

plot(((AT$V3-AT$V2)/(tunit$V3-tunit$V2)), log10(AT$V3-AT$V2))
sum((AT$V3-AT$V2)/(tunit$V3-tunit$V2) <= 0.2) / length(AT$V3)


## use AT that <=20% of tunit length
navg=100
heatmap.AT4<-function(at, t, r=0.2){
tunit <-  read.table(paste("../",at,"_AT_3tunitIntersectNativeHMM_tunits.bed", sep=""), header = F)
AT <- read.table(paste("../",at,"_AT_3tunitIntersectNativeHMM_intersectRegion.bed", sep=""), header = F)
AT$AT.tunit.length.ratio <- (AT$V3-AT$V2)/(tunit$V3-tunit$V2)
#file.plus.bw <- paste("../",at,"_all_plus.rpm.bw", sep="")
#file.minus.bw <- paste("../",at,"_all_minus.rpm.bw", sep="")
  file.bw.plus.pat <- paste("../",t,"_MB6_all_R1.pat_1bp_plus.bw", sep="")
  file.bw.minus.pat <- paste("../",t,"_MB6_all_R1.pat_1bp_minus.bw", sep="")
  file.bw.plus.mat <- paste("../",t,"_MB6_all_R1.mat_1bp_plus.bw", sep="")
  file.bw.minus.mat <- paste("../",t,"_MB6_all_R1.mat_1bp_minus.bw", sep="")
  hl.bw.plus.pat <- paste("../",at,"_MB6_all_R1.pat_1bp_plus.bw", sep="")
  hl.bw.minus.pat <- paste("../",at,"_MB6_all_R1.pat_1bp_minus.bw", sep="")
  hl.bw.plus.mat <- paste("../",at,"_MB6_all_R1.mat_1bp_plus.bw", sep="")
  hl.bw.minus.mat <- paste("../",at,"_MB6_all_R1.mat_1bp_minus.bw", sep="")

#AT <- AT[AT$AT.tunit.ratio <= 0.2,]
heatmap.AT3(AT[AT$AT.tunit.length.ratio <= r,], file.bw.plus.pat,file.bw.minus.pat, file.bw.plus.mat ,file.bw.minus.mat, 
            hl.bw.plus.pat=hl.bw.plus.pat ,hl.bw.minus.pat=hl.bw.minus.pat, hl.bw.plus.mat = hl.bw.plus.mat ,hl.bw.minus.mat=hl.bw.minus.mat, 
            dist=20000, step=500,  
            file.pdf=paste("twentyPercent_",at,"-AT_",t,"-bw_break5",".pdf",sep=""), 
            bl_wd=1, show.AT.line=TRUE, navg=navg, use.log=FALSE, times=1, 
            breaks = seq(0,5,0.01))
heatmap.AT3(AT[AT$AT.tunit.length.ratio <= r,], file.bw.plus.pat,file.bw.minus.pat, file.bw.plus.mat ,file.bw.minus.mat, 
            hl.bw.plus.pat=hl.bw.plus.pat ,hl.bw.minus.pat=hl.bw.minus.pat, hl.bw.plus.mat = hl.bw.plus.mat ,hl.bw.minus.mat=hl.bw.minus.mat, 
            dist=20000, step=500,  
            file.pdf=paste("twentyPercent_",at,"-AT_",t,"-bw_break2",".pdf",sep=""), 
            bl_wd=1, show.AT.line=TRUE, navg=navg, use.log=FALSE, times=1, 
            breaks = seq(0,2,0.01))
}
for (at in c("BN","LV")){
  for (t in c("BN","LV")){
    heatmap.AT4(at, t, r=0.2)
  }
}

