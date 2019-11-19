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

effect="strain+imprinting.bed"
for (effect in c("strain+imprinting.bed","imprinting.bed","strain.bed")){
for (OR in c("BN", "SP", "HT", "SK", "KD", "ST", "GI", "LV")){
  pdf(paste(OR, "TSSinCluster",effect,"seperate.pdf", sep = "_"))
  par(mar=c(6.1, 7.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
  par(mgp=c(3,1,0))
  par(cex.lab=2.2, cex.axis=2.2)
#OR="LV"
f1=read.table(paste(OR,"dREG_count_in_T8_2Strand_p0.05_effect", effect, "cluster", sep = "_"),header=F)
f2=read.table(paste(OR,"dREG_count_in_T8_2Strand_p0.05_effect", effect, "cluster_with", sep = "_"),header=F)
f3=read.table(paste(OR,"dREG_count_in_T8_2Strand_p0.05_effect", effect, "cluster_without", sep = "_"),header=F)


hist(f3$V1,col="red", density=25
     , breaks = seq(0,200,binSize)
     , freq = F
     ,ylab="Proportion of domains"
     , xlim=c(0,50)
     ,xlab=paste("Number of", OR,"TSSs in each domain",sep=" ")
     ,main= ""
     ,add=F
     ,las=1
)


hist(f2$V1,col="blue"
     , breaks = seq(0,200,binSize)
     , freq = F
     #, xlab="HMM blocks in the domain"
     #, add=T
     , xlim=c(0,50)
     ,xlab="Number of TSSs in each domain"
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
pdf(paste(OR, "TSSinCluster", effect, ".pdf", sep = "_"))
par(mar=c(6.1, 7.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
par(mgp=c(3,1,0))
par(cex.lab=2.2, cex.axis=2.2)
hist(f1$V1,col="dark green", density=25
     , breaks = seq(0,200,binSize)
     , freq = F
     ,ylab="Proportion of domains"
     , xlim=c(0,50)
     ,xlab=paste("Number of", OR,"TSSs in each domain",sep=" ")
     ,main= ""
     ,add=F
     ,las=1
)
dev.off()
}
}

for (OR in c("BN", "SP", "HT", "SK", "KD", "ST", "GI", "LV")){
  pdf(paste(OR, "TSSinCluster_seperate.pdf", sep = "_"))
  par(mar=c(6.1, 7.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
  par(mgp=c(3,1,0))
  par(cex.lab=2.2, cex.axis=2.2)
binSize=1

f2=read.table(paste(OR,"_dREG_count_in_T8_2Strand_p0.05_effect_strain.bed_cluster", sep = ""),header=F)
a=hist(f2$V1,col="blue" 
       #,density=25, 
       ,breaks = seq(0,200,binSize)
       , freq = F
       ,ylab="Proportion of domains"
       , xlim=c(0,50)
       ,xlab=paste("Number of","TSSs in each domain",sep=" ")
       ,main= ""
       ,add=F
       ,las=1
)
f1=read.table(paste(OR, "_dREG_count_in_T8_2Strand_p0.05_effect_imprinting.bed_cluster", sep = ""),header=F)
a=hist(f1$V1,col="red" 
       ,density=25
       , breaks = seq(0,200,binSize)
       , freq = F
       ,ylab="Proportion of domains"
       , xlim=c(0,50)
       ,xlab=paste("Number of","TSSs in each domain",sep=" ")
       ,main= ""
       ,add=T
       ,las=1
)
legend("topright", 
       legend = c( "Imprinting","Strain effect"), 
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
      
      a=hist(f1$V1,col="dark green", density=25
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
binSize=1
f=read.table("gencode.vM20.annotation_geneMerged.bed_count_in_T8_2Strand_p0.05_effect_strain+imprinting.bed_cluster",header=F)
a=hist(f$V1,col="dark green"#, density=25
     , breaks = seq(0,200,binSize)
     , freq = F
     ,ylab="Proportion of clusters"
     , xlim=c(0,50)
     ,xlab=paste("Number of","gencode gene annotations in each cluster",sep=" ")
     ,main= ""
     ,add=F
     ,las=1
)
plot(f$V1, f$V4-f$V3, xlab = "Number of gencode gene annotations in each cluster", ylab = "cluster length", las=1)

hist(f$V1[f$V4-f$V3>1000000],col="red", density=25
     , breaks = seq(0,200,binSize)
     , freq = F
     ,ylab="Proportion of clusters"
     , xlim=c(0,50)
     ,xlab=paste("Number of","gencode gene annotations in each cluster",sep=" ")
     ,main= ""
     ,add=T
     ,las=1
)
hist(f$V1[f$V4-f$V3>100000],col="blue", density=25
     , breaks = seq(0,200,binSize)
     , freq = F
     ,ylab="Proportion of clusters"
     , xlim=c(0,50)
     ,xlab=paste("Number of","gencode gene annotations in each cluster",sep=" ")
     ,main= ""
     ,add=T
     ,las=1
)
legend("topright", 
       legend = c( "all cluster","cluster > 100Kb","cluster > 1Mb"), 
       #pch=c(15,15),
       cex=2, 
       lty=c(0,0),
       #bty="n",
       lwd=1.5, 
       density=c(10000, 25, 25),
       angle=c(180,45,45),
       #angle=45,
       fill=c("dark green","blue","red")
       , bty = "n"
)

dev.off()

pdf("gencode.vM20.annotation_geneMergedinCluster_SI.pdf")
par(mar=c(6.1, 7.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
par(mgp=c(3,1,0))
par(cex.lab=2.2, cex.axis=2.2)
binSize=1
f2=read.table("gencode.vM20.annotation_geneMerged.bed_count_in_T8_2Strand_p0.05_effect_strain.bed_cluster",header=F)
f2_0=read.table("T8_2Strand_p0.05_effect_strain.bed_cluster",header=F)

a=hist(c(rep(0,dim(f2_0)[1] - dim(f2)[1]),f2$V1),col="blue" 
       #,density=25, 
       ,breaks = seq(0,200,binSize)
       , freq = F
       ,ylab="Proportion of domains"
       , xlim=c(0,50)
       ,xlab=paste("Number of","gencode gene annotations in each domain",sep=" ")
       ,main= ""
       ,add=F
       ,las=1
)
f1=read.table("gencode.vM20.annotation_geneMerged.bed_count_in_T8_2Strand_p0.05_effect_imprinting.bed_cluster",header=F)
f1_0=read.table("T8_2Strand_p0.05_effect_imprinting.bed_cluster",header=F)

a=hist(c(rep(0,dim(f1_0)[1] - dim(f1)[1]),f1$V1),col="red" 
       ,density=25
       , breaks = seq(0,200,binSize)
       , freq = F
       ,ylab="Proportion of clusters"
       , xlim=c(0,50)
       ,xlab=paste("Number of","gencode gene annotations in each cluster",sep=" ")
       ,main= ""
       ,add=T
       ,las=1
)
legend("topright", 
       legend = c( "Imprinting","Strain effect"), 
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
plot(a$mids,log10(a$counts))


setwd("~/Box Sync/Danko_lab_work/F1_8Tissues/Cluster")

f=read.table("T8_2Strand_p0.05_effect_strain+imprinting.bed_cluster_length",header=F)
hist(log10(f$V1),col="dark green", density=25
     #, breaks = seq(0,10000000,100000)
     , freq = F
     ,ylab="Number of clusters"
     #, xlim=c(0,50)
     ,xlab="Cluster length"
     ,main= ""
     ,add=F
     ,las=2
)
f2=read.table("T8_2Strand_p0.05_effect_strain.bed_cluster_length",header=F)
hist(log10(f2$V1),col="blue" 
     #,density=25
     , breaks = seq(0,8,1)
     , freq = F
     ,ylab="Percentage of domains"
     , xlim=c(0,7)
     ,xlab="Domain length"
     ,main= ""
     ,add=T
     ,las=2
)

f1=read.table("T8_2Strand_p0.05_effect_imprinting.bed_cluster_length",header=F)
hist(log10(f1$V1),col="red" 
     ,density=25
     , breaks = seq(0,8,1)
     , freq = F
     ,ylab="Percentage of domains"
     , xlim=c(0,7)
     ,xlab="Domain length"
     ,main= ""
     ,add=T
     ,las=2
)
legend("topleft", 
       legend = c( "Imprinting","Strain effect"), 
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


