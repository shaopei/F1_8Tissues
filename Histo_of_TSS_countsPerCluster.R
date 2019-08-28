setwd("~/Box Sync/Danko_lab_work/F1_8Tissues/Cluster/TSSinCluster/")

###color transparant
t_col <- function(color, percent = 50, name = NULL) {
  #      color = color name
  #       percent = % transparency
  #          name = an optional name for the color
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100-percent)*255/100,
               names = name)
  ## Save the color
  invisible(t.col)
}

mycol1 <- t_col("blue", perc = 50, name = "lt.blue")
mycol2 <- t_col("green", perc = 50, name = "lt.green")
binSize=1

for (OR in c("BN", "SP", "HT", "SK", "KD", "ST", "GI", "LV")){
  pdf(paste(OR, "TSSinCluster_seperate.pdf", sep = "_"))
  par(mar=c(6.1, 7.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
  par(mgp=c(3,1,0))
  par(cex.lab=2.2, cex.axis=2.2)
#OR="LV"
f1=read.table(paste(OR,"dREG_count_in_T8_2Strand_p0.05_effect_strain+imprinting.bed_cluster", sep = "_"),header=F)
f2=read.table(paste(OR,"dREG_count_in_T8_2Strand_p0.05_effect_strain+imprinting.bed_cluster_with", sep = "_"),header=F)
f3=read.table(paste(OR,"dREG_count_in_T8_2Strand_p0.05_effect_strain+imprinting.bed_cluster_without", sep = "_"),header=F)

hist(f3$V1,col="red", density=25
     , breaks = seq(0,200,binSize)
     , freq = F
     ,ylab="Proportion of clusters"
     , xlim=c(0,50)
     ,xlab=paste("Number of", OR,"TSSs in each cluster",sep=" ")
     ,main= ""
     ,add=F
     ,las=1
)


hist(f2$V1,col="blue"
     , breaks = seq(0,200,binSize)
     , freq = F
     #, xlab="HMM blocks in the cluster"
     #, add=T
     , xlim=c(0,50)
     ,xlab="Number of TSSs in each cluster"
     ,main=""
     ,add=T
)
hist(f3$V1,col="red", density=25
     , breaks = seq(0,200,binSize)
     , freq = F
     ,add=T
)


legend("topright", 
       legend = c( paste(OR,"Allelic Balance",sep=" "),  paste(OR,"Allelic Biased",sep=" ")),
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
pdf(paste(OR, "TSSinCluster.pdf", sep = "_"))
par(mar=c(6.1, 7.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
par(mgp=c(3,1,0))
par(cex.lab=2.2, cex.axis=2.2)
hist(f1$V1,col="dark green", density=25
     , breaks = seq(0,200,binSize)
     , freq = F
     ,ylab="Proportion of clusters"
     , xlim=c(0,50)
     ,xlab=paste("Number of", OR,"TSSs in each cluster",sep=" ")
     ,main= ""
     ,add=F
     ,las=1
)
dev.off()
}

setwd("~/Box Sync/Danko_lab_work/F1_8Tissues/Cluster/AlleleHMMInCluster/")

for (OR in c("BN", "SP", "HT", "SK", "KD", "ST", "GI", "LV")){
  for (cross in c("MB6", "PB6")){
    for (s in c("plus", "minus")){
      pdf(paste("pdf/",paste(OR,cross,"HMM",s, "HMMInCluster.pdf", sep = "_"), sep=""))
      par(mar=c(6.1, 7.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
      par(mgp=c(3,1,0))
      par(cex.lab=2.2, cex.axis=2.2)
      f1=read.table(paste(OR,cross,"HMM",s,"count_in_T8_2Strand_p0.05_effect_strain+imprinting.bed_cluster", sep = "_"),header=F)
      
      hist(f1$V1,col="dark green", density=25
           , breaks = seq(0.5,max(10.1,max(f1$V1)+0.5),1)
           #, freq = F
           ,ylab="Number of clusters"
           #, xlim=c(0,50)
           ,xlab=paste("Number of", OR,s,"AlleleHMM blocks in each cluster",sep=" ")
           ,main= ""
           ,add=F
           ,las=1
      )
      dev.off()
    }
  }
}

  
  
  
setwd("~/Box Sync/Danko_lab_work/F1_8Tissues/Cluster/GeneAnnotationInCluster/")
pdf("gencode.vM20.annotation_geneMergedinCluster.pdf")
par(mar=c(6.1, 7.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
par(mgp=c(3,1,0))
par(cex.lab=2.2, cex.axis=2.2)
f=read.table("gencode.vM20.annotation_geneMerged.bed_count_in_T8_2Strand_p0.05_effect_strain+imprinting.bed_cluster",header=F)
a=hist(f$V1,col="dark green", density=25
     , breaks = seq(0,200,binSize)
     , freq = F
     ,ylab="Proportion of clusters"
     , xlim=c(0,50)
     ,xlab=paste("Number of","gencode gene annotations in each cluster",sep=" ")
     ,main= ""
     ,add=F
     ,las=1
)
dev.off()
plot(a$mids,log10(a$counts))



f=read.table("T8_2Strand_p0.05_effect_strain+imprinting.bed_cluster_length",header=F)
hist(log10(f$V1),col="dark green", density=25
     #, breaks = seq(0,10000000,100000)
     #, freq = F
     ,ylab="Number of clusters"
     #, xlim=c(0,50)
     ,xlab="Cluster length (log10)"
     ,main= ""
     ,add=F
     ,las=1
)



##END

lamda = with(f2, 1/mean(V1))
x=seq(0,200,1)
lines(x, lamda*exp(-lamda*x), col="dark blue",lwd=1.5)
lamda = with(f3, 1/mean(V1))
x=seq(0,200,1)
lines(x, lamda*exp(-lamda*x), col="dark red",lwd=1.5)

#lines(density(f2$V1), col="light blue") 
legend("topright", 
       legend = paste("Exponential pdf lamda = ", format(lamda, digits=2, nsmall=2),sep=""),
       lty=1, lwd=1.5, col="orange", bty = "n")

## END


