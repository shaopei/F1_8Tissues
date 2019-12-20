setwd("~/Box Sync/Danko_lab_work/F1_8Tissues/PolyA_Allele-specific/heatmap/")


read_count <- function(file.plus.bw, file.minus.bw , bed6) {
    bw.plus  <- load.bigWig( file.plus.bw )
    bw.minus <- load.bigWig( file.minus.bw )
    CountMatrix <- bed6.region.bpQuery.bigWig(bw.plus, bw.minus, bed6[,c(1:6)] ,abs.value=TRUE, op = "sum")
    unload.bigWig(bw.plus);
    unload.bigWig(bw.minus);
    
    return(CountMatrix);    
}

Tunit_BW_pairs <- function (at, t){
tunit <-  read.table(paste("../",at,"_AT_3tunitIntersectNativeHMM_tunits.bed", sep=""), header = F)
file.plus.bw <- paste("../",t,"_all_plus.rpm.bw", sep="")
file.minus.bw <- paste("../",t,"_all_minus.rpm.bw", sep="")
log(read_count(file.plus.bw, file.minus.bw, tunit)+1)
}

BN.Tunit_BN.bw <- Tunit_BW_pairs ("BN","BN")
BN.Tunit_LV.bw <- Tunit_BW_pairs ("BN","LV")
LV.Tunit_LV.bw <- Tunit_BW_pairs ("LV","LV")
LV.Tunit_BN.bw <- Tunit_BW_pairs ("LV","BN")
summary(BN.Tunit_BN.bw)

library("vioplot")
library(RColorBrewer)
display.brewer.all()
display.brewer.pal(n = 9,name = "PuOr")
PuOr = brewer.pal(n = 9,name = "PuOr")


BN_col = PuOr[7]
LV_col = PuOr[3] #'#E69F00'
BN_l = PuOr[8]
LV_l = PuOr[2] #'#E69F00'
vioplot(BN.Tunit_BN.bw, BN.Tunit_LV.bw, LV.Tunit_BN.bw, LV.Tunit_LV.bw,
        names=c("BN.BN", "BN.LV", "LV.BN", "LV.LV"), 
        border=c(BN_l, BN_l, LV_l, LV_l), lwd=4,
        col=c(BN_col, LV_col,BN_col, LV_col),
        xlab="Tunit.BW",
        ylab="log(expression)",
        ylim = c(-1,9), las=2)

vioplot( LV.Tunit_BN.bw, LV.Tunit_LV.bw,BN.Tunit_BN.bw, BN.Tunit_LV.bw,
        names=c( "LV.BN", "LV.LV", "BN.BN", "BN.LV"),
        border=c(LV_l, LV_l,BN_l, BN_l), lwd=4,
        col=c(BN_col, LV_col,BN_col, LV_col),
        xlab="Tunit.BW",
        ylab="log(expression)",
        ylim = c(-1,9), las=2)

legend("topleft", legend=c("Brain", "Liver"),
       fill=c(BN_col, LV_col), cex=1.5 , bty = "n",
       border = c(BN_l, LV_l))

