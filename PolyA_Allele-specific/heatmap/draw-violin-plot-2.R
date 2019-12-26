#domain_length.pdf
#setwd("~/Box Sync/Danko_lab_work/F1_8Tissues/PolyA_Allele-specific/")

args=(commandArgs(TRUE))

folder = args[1]
f1_p = args[2] #"BN_AT_3tunitIntersectNativeHMM_intersectRegion.bed"
f2_p = args[3] #"BN_AT_3tunitIntersectNativeHMM_tunit.bed"
t=args[4]

setwd(folder) #("/Volumes/SPC_SD/IGV/PolyA_Allele-specific")
pdf(paste(t, "length.pdf", sep = "_"))
par(mar=c(6.1, 7.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
par(mgp=c(3,1,0))
par(cex.lab=2.2, cex.axis=2.2)
f1=read.table(f1_p,header=F)
h<- hist(log10(f1$V3-f1$V2),#col="red" 
#         ,density=25
         , breaks = seq(0,7,1)
#         , freq = T
#         , prob=TRUE
#         ,ylab="Proportion of domains"
#         , xlim=c(0,7)
#         ,xlab="AlleleHMM block length (log10)"
#         ,main= t
#         ,add=F
#         ,las=2
         ,plot =FALSE
         ,right = FALSE
)

h$counts=h$counts/sum(h$counts)

f2=read.table(f2_p,header=F)
h2<- hist(log10(f2$V3-f2$V2),col="blue" 
          , breaks = seq(0,8,1)
          , freq = F
          #,add=F
          #,las=2
          ,plot =FALSE
          ,right = FALSE
)

h2$counts=h2$counts/sum(h2$counts)
plot(h,col="red" 
     ,density=25     
     ,ylab="Proportion"
     , xlim=c(0,7)
     ,ylim=c(0,max(h$density,h2$density))
     ,xlab="block length"     
     ,las=2
     ,xaxt='n',main= t)
axis(1, at=seq(0,7,1), labels=c(0,10,100,1000,"10,000","100,000","1000,000","10,000,000"), las=2)

plot(h2,col="blue" , add=T)
plot(h,col="red" ,density=25, add=T)

legend("topleft", 
       legend = c( "AT window","Tunit"), 
       #pch=c(15,15),
       cex=2, 
       lty=c(0,0),
       #bty="n",
       lwd=1.5, 
       density=c(25, 10000),
       angle=c(45, 180),
       #angle=45,
       fill=c("red","blue")
       , bty = "n"
)
dev.off()
####

library("vioplot")
library(bigWig)

read_count <- function(file.plus.bw, file.minus.bw , bed6) {
  bw.plus  <- load.bigWig( file.plus.bw )
  bw.minus <- load.bigWig( file.minus.bw )
  CountMatrix <- bed6.region.bpQuery.bigWig(bw.plus, bw.minus, bed6[,c(1:6)] ,abs.value=TRUE, op = "sum")
  unload.bigWig(bw.plus);
  unload.bigWig(bw.minus);
  
  return(CountMatrix);    
}

#tussue_list <-c("BN", "HT", "SK", "SP", "KD", "LV", "GI", "ST")
tussue_list <-c("BN",  "SP", "LV")


par(mfrow=c(1,2))
Tversion <- "4tunitIntersectNativeHMM"
for (Tversion in c("3tunitIntersectNativeHMM","4tunitIntersectNativeHMM")){
  
getTunitLength <- function(t, path="./", body=paste("_AT",Tversion, "tunits.bed", sep = "_")){
  df <- read.table(paste(path,t,body, sep=""), header = F)
  return(log10(df$V3-df$V2))
}

getATLength <- function(t, path="./", body=paste("_AT",Tversion, "intersectRegion.bed", sep = "_")){
  df <- read.table(paste(path,t,body, sep=""), header = F)
  return(log10(df$V3-df$V2))
} 

setwd(paste("~/Box Sync/Danko_lab_work/F1_8Tissues/PolyA_Allele-specific/tunitIntersectNativeHMM",Tversion,sep = "/"))
vioplot(getTunitLength(tussue_list[1]),getTunitLength(tussue_list[2]),getTunitLength(tussue_list[3]), #getTunitLength(tussue_list[4]), getTunitLength(tussue_list[5]),getTunitLength(tussue_list[6]),getTunitLength(tussue_list[7]),getTunitLength(tussue_list[8]),
        names = tussue_list ,
        ylab= "log10(Tunit length)", main=Tversion )

vioplot(getATLength(tussue_list[1]),getATLength(tussue_list[2]),getATLength(tussue_list[3]),#getATLength(tussue_list[4]),getATLength(tussue_list[5]),getATLength(tussue_list[6]),getATLength(tussue_list[7]),getATLength(tussue_list[8]),
        names = tussue_list ,
        ylab= "log10(AT length)",
        ylim = c(0,7),
        main=Tversion )
}

### the relationship between Tunit expression level and AT window length
par(mfrow=c(2,3))
for (at in c("BN", "SK", "SP", "LV", "GI", "ST")){
#for (at in c("BN",  "SP", "LV")){
tunit <-  read.table(paste("./",at,"_AT_",Tversion, "_tunits.bed", sep=""), header = F)
AT <- read.table(paste("./",at,"_AT_",Tversion,"_intersectRegion.bed", sep=""), header = F)
file.plus.bw <- paste("/Volumes/SPC_SD/IGV/bigWig_notBN/",at,"_all_plus.rpm.bw", sep="")
file.minus.bw <- paste("/Volumes/SPC_SD/IGV/bigWig_notBN/",at,"_all_minus.rpm.bw", sep="")
tunit_rpkm <- read_count(file.plus.bw, file.minus.bw, tunit) *1000/(tunit$V3 - tunit$V2) 
AT_rpkm <- read_count(file.plus.bw, file.minus.bw, AT) *1000/(AT$V3 - AT$V2) 
plot(log10(tunit_rpkm), log10(AT_rpkm), main=paste(at,Tversion,sep=" "), xlim = c(-2,2), ylim = c(-2,2))
#tunit_rpmSum <- read_count(file.plus.bw, file.minus.bw, tunit) 
#plot_colorByDensity(log10(AT$V3-AT$V2), log(tunit_rpkm), main=paste(at,Tversion,sep=" "),
#                    xlab="log10(AT$V3-AT$V2)", ylab=" log(tunit_rpkm)")
#smoothScatter(log10(AT$V3-AT$V2), log(tunit_rpkm), main=paste(at,Tversion,sep=" "))
#points(log10(AT$V3-AT$V2), log(tunit_rpkm))
# plot_colorByDensity(log10(AT$V3-AT$V2), log(tunit_rpmSum), main=paste(at,Tversion,sep=" "),
#                    xlab="log10(AT$V3-AT$V2)", ylab=" log(tunit_rpmSum)")
#smoothScatter(log(AT$V3-AT$V2), log(tunit_rpmSum), main=paste(at,Tversion,sep=" "))
#points(log(AT$V3-AT$V2), log(tunit_rpmSum))
}

# 
plot_colorByDensity = function(x1,x2,
                               ylim=c(min(x2),max(x2)),
                               xlim=c(min(x1),max(x1)),
                               xlab="",ylab="",main="") {

  df <- data.frame(x1,x2)
  x <- densCols(x1,x2, colramp=colorRampPalette(c("black", "white")))
  df$dens <- col2rgb(x)[1,] + 1L
  cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F","#FCFF00", "#FF9400", "#FF3100"))(256)
  df$col <- cols[df$dens]
  plot(x2~x1, data=df[order(df$dens),],
       ylim=ylim,xlim=xlim,pch=20,col=col,
       cex=2,xlab=xlab,ylab=ylab,
       main=main)
}

### plots of expression level

Tunit_BW_pairs <- function (at, t){
  tunit <-  read.table(paste("./",at,"_AT_",Tversion, "_tunits.bed", sep=""), header = F)
  file.plus.bw <- paste("/Volumes/SPC_SD/IGV/bigWig_notBN/",t,"_all_plus.rpm.bw", sep="")
  file.minus.bw <- paste("/Volumes/SPC_SD/IGV/bigWig_notBN/",t,"_all_minus.rpm.bw", sep="")
  log10(read_count(file.plus.bw, file.minus.bw, tunit)*1000/(tunit$V3 - tunit$V2)+0.01) # log10(rpkm+0.01)
}
AT_BW_pairs <- function (at, t){
  AT <- read.table(paste("./",at,"_AT_",Tversion,"_intersectRegion.bed", sep=""), header = F)
  file.plus.bw <- paste("/Volumes/SPC_SD/IGV/bigWig_notBN/",t,"_all_plus.rpm.bw", sep="")
  file.minus.bw <- paste("/Volumes/SPC_SD/IGV/bigWig_notBN/",t,"_all_minus.rpm.bw", sep="")
  log10(read_count(file.plus.bw, file.minus.bw, AT)*1000/(AT$V3 - AT$V2)+0.01) # log10(rpkm+0.01)
}



dev.off()
par(mfrow=c(1,9))
for (at in c("BN",  "SP", "LV")){
  for (t in c("BN",  "SP", "LV")){
    #vioplot(Tunit_BW_pairs (at,t), ylim = c(-2,2),names=paste(at,".Tu_",t,".bw", sep=""))
    vioplot(AT_BW_pairs (at,t), ylim = c(-2,2),names=paste(at,".AT_",t,".bw", sep=""))
  }
}

Bed_BW_pairs <- function (at, t, path="./", body=paste("_AT_",Tversion,"_intersectRegion.bed", sep="")){
  Bed <- read.table(paste(path,at,body, sep=""), header = F)
  file.plus.bw <- paste("/Volumes/SPC_SD/IGV/bigWig_notBN/",t,"_all_plus.rpm.bw", sep="")
  file.minus.bw <- paste("/Volumes/SPC_SD/IGV/bigWig_notBN/",t,"_all_minus.rpm.bw", sep="")
  log10(read_count(file.plus.bw, file.minus.bw, Bed)*1000/(Bed$V3 - Bed$V2)+0.001) # log10(rpkm+0.01)
}

par(mfrow=c(3,3))
for (t in c("BN",  "SP", "LV")){
  for (at in c("BN",  "SP", "LV")){
    if (t!=at){
    col="white"
    }else{
      col="red"
    }
    x <- Bed_BW_pairs (at,t, body=paste("_AT_",Tversion,"_intersectRegion.bed", sep=""))
    hist(x,
         main = paste(at,".AT_",t,".bw", spe=""),
         breaks = seq(min(x)-0.1, max(x)+0.1,0.1),
         ylim=c(0,100), xlim = c(-3,3), col=col 
         
         )#, breaks = seq(-2,2,0.1))#, border = BN_l, col = BN_col, ylim = c(0,100))
  }
}
for (t in c("BN",  "SP", "LV")){
  for (at in c("BN",  "SP", "LV")){
    if (t!=at){
      col="white"
    }else{
      col="blue"
    }
    hist(Bed_BW_pairs (at,t,body=paste("_AT_",Tversion,"_tunits.bed", sep="")),
         main = paste(at,".Tunit_",t,".bw", spe=""),
         breaks = seq(-3,3,0.1),
         ylim=c(0,100), xlim = c(-3,3), col=col 
         
    )#, breaks = seq(-2,2,0.1))#, border = BN_l, col = BN_col, ylim = c(0,100))
  }
}

### absolute expression level vs relative level of native organ and other organ
par(mfrow=c(3,3))
for (t in c("BN",  "SP", "LV")){
  for (at in c("BN",  "SP", "LV")){
    if (t!=at){
      col="black"
    }else{
      col="red"
    }
    plot(AT_BW_pairs(at,at),AT_BW_pairs(at,t) ,
         ylab=paste(at,".AT_",t,".bw", spe=""), xlab=paste(at,".AT_",at,".bw", spe=""),
         col=col)
  }
}
for (t in c("BN",  "SP", "LV")){
  for (at in c("BN",  "SP", "LV")){
    if (t!=at){
      col="black"
    }else{
      col="blue"
    }
    plot(Tunit_BW_pairs(at,at),Tunit_BW_pairs(at,t) ,
         ylab=paste(at,".Tunit_",t,".bw", spe=""), xlab=paste(at,".Tunit_",at,".bw", spe=""),
         col=col)
  }
}


# absolute expression 
navg=1
#at="LV"; other = "BN"
#at="BN"; t = "LV"

tissues <- c("BN", "LV","SP")
for (at in tissues){
  at="SP"
  for (t in tissues){
    AT <- read.table(paste("./",at,"_AT_",Tversion,"_intersectRegion.bed", sep=""), header = F)
    t.value <-Tunit_BW_pairs (at, t)
    
    
    #for (t in c("LV", "BN")){
    file.bw.plus.pat <- paste("/Volumes/SPC_SD/IGV/bigWig_notBN/",t,"_MB6_all_R1.pat_1bp_plus.bw", sep="")
    file.bw.minus.pat <- paste("/Volumes/SPC_SD/IGV/bigWig_notBN/",t,"_MB6_all_R1.pat_1bp_minus.bw", sep="")
    file.bw.plus.mat <- paste("/Volumes/SPC_SD/IGV/bigWig_notBN/",t,"_MB6_all_R1.mat_1bp_plus.bw", sep="")
    file.bw.minus.mat <- paste("/Volumes/SPC_SD/IGV/bigWig_notBN/",t,"_MB6_all_R1.mat_1bp_minus.bw", sep="")
    hl.bw.plus.pat <- paste("/Volumes/SPC_SD/IGV/bigWig_notBN/",at,"_MB6_all_R1.pat_1bp_plus.bw", sep="")
    hl.bw.minus.pat <- paste("/Volumes/SPC_SD/IGV/bigWig_notBN/",at,"_MB6_all_R1.pat_1bp_minus.bw", sep="")
    hl.bw.plus.mat <- paste("/Volumes/SPC_SD/IGV/bigWig_notBN/",at,"_MB6_all_R1.mat_1bp_plus.bw", sep="")
    hl.bw.minus.mat <- paste("/Volumes/SPC_SD/IGV/bigWig_notBN/",at,"_MB6_all_R1.mat_1bp_minus.bw", sep="")
    file.plus.bw <- paste("/Volumes/SPC_SD/IGV/bigWig_notBN/",t,"_all_plus.rpm.bw", sep="")
    file.minus.bw <- paste("/Volumes/SPC_SD/IGV/bigWig_notBN/",t,"_all_minus.rpm.bw", sep="")
    
    heatmap.AT3(AT[t.value > -1,], file.bw.plus.pat,file.bw.minus.pat, file.bw.plus.mat ,file.bw.minus.mat,
                hl.bw.plus.pat=hl.bw.plus.pat ,hl.bw.minus.pat=hl.bw.minus.pat, hl.bw.plus.mat = hl.bw.plus.mat ,hl.bw.minus.mat=hl.bw.minus.mat,
                dist=20000, step=500,
                file.pdf=paste(at,".AT_",t,".bw_",1,".pdf",sep=""),
                bl_wd=1, show.AT.line=TRUE, navg=navg, use.log=FALSE,
                times=10, breaks = seq(0,5,0.01))
    heatmap.AT3(AT[t.value <= -1,], file.bw.plus.pat,file.bw.minus.pat, file.bw.plus.mat ,file.bw.minus.mat,
                hl.bw.plus.pat=hl.bw.plus.pat ,hl.bw.minus.pat=hl.bw.minus.pat, hl.bw.plus.mat = hl.bw.plus.mat ,hl.bw.minus.mat=hl.bw.minus.mat,
                dist=20000, step=500,
                file.pdf=paste(at,".AT_",t,".bw_",2,".pdf",sep=""),
                bl_wd=1, show.AT.line=TRUE, navg=navg, use.log=FALSE,
                times=1000, breaks = seq(0,5,0.01))
    # heatmap.AT4_allreads(AT[t.value > -1,], file.plus.bw,file.minus.bw, file.plus.bw ,file.minus.bw, 
    #                      dist=20000, step=500,
    #                      file.pdf=paste(t,".bw_",1,".pdf",sep=""), 
    #                      bl_wd=1, show.AT.line=TRUE, navg=navg, times=30, use.log=FALSE, breaks=seq(0, 10, 1))
    # heatmap.AT4_allreads(AT[t.value <= -1,], file.plus.bw,file.minus.bw, file.plus.bw ,file.minus.bw, 
    #                      dist=20000, step=500,
    #                      file.pdf=paste(t,".bw_",2,".pdf",sep=""), 
    #                      bl_wd=1, show.AT.line=TRUE, navg=navg, times=30, use.log=FALSE, breaks=seq(0, 10, 1))
     
  }
}

