#domain_length.pdf
setwd("~/Box Sync/Danko_lab_work/F1_8Tissues/Cluster")
pdf("domain_length.pdf")
par(mar=c(6.1, 7.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
par(mgp=c(3,1,0))
par(cex.lab=2.2, cex.axis=2.2)
f1=read.table("T8_2Strand_p0.05_effect_imprinting.bed_cluster_length",header=F)
h<- hist(log10(f1$V1),col="red" 
     ,density=25
     , breaks = seq(0,8,0.5)
     #, freq = F
     , prob=TRUE
     ,ylab="Proportion of domains"
     , xlim=c(0,7)
     ,xlab="Domain length"
     ,main= ""
     ,add=F
     ,las=2
     ,plot =FALSE
)

h$counts=h$counts/sum(h$counts)
plot(h,col="red" 
     ,density=25     
     ,ylab="Proportion of domains"
     , xlim=c(0,7)
     ,xlab="Domain length"     
     ,las=2
     ,xaxt='n',main= "")
axis(1, at=seq(0,7,1), labels=c(0,10,100,1000,"10,000","100,000","1000,000","10,000,000"), las=2)


f2=read.table("T8_2Strand_p0.05_effect_strain.bed_cluster_length",header=F)
h2<- hist(log10(f2$V1),col="blue" 
     , breaks = seq(0,8,0.5)
     , freq = F
     #,add=F
     #,las=2
     ,plot =FALSE
)
h2$counts=h2$counts/sum(h2$counts)
plot(h2,col="blue" , add=T)
plot(h,col="red" ,density=25, add=T)

legend("topleft", 
       legend = c( "Imprinted","Strain effect"), 
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




# gencode.vM20.annotation_geneMergedinCluster_SI.pdf
setwd("~/Box Sync/Danko_lab_work/F1_8Tissues/Cluster/GeneAnnotationInCluster/")
pdf("gencode.vM20.annotation_geneMergedinCluster_SI.pdf")
par(mar=c(6.1, 7.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
par(mgp=c(3,1,0))
par(cex.lab=2.2, cex.axis=2.2)
binSize=1
f2=read.table("gencode.vM20.annotation_geneMerged.bed_count_in_T8_2Strand_p0.05_effect_strain.bed_cluster",header=F)
f2_0=read.table("T8_2Strand_p0.05_effect_strain.bed_cluster",header=F)

h2=hist(c(rep(0,dim(f2_0)[1] - dim(f2)[1]),f2$V1),col="blue" 
       #,density=25, 
       ,breaks = seq(0,200,binSize)
       , freq = F
       ,ylab="Proportion of domains"
       , xlim=c(0,50)
       ,xlab=paste("Number of","gencode gene annotations in each domain",sep=" ")
       ,main= ""
       ,add=F
       ,las=1
       ,plot =FALSE
)

h2$counts=h2$counts/sum(h2$counts)
plot(h2,col="blue" 
     ,ylab="Proportion of domains"
     ,xlab=paste("Number of","gencode gene annotations in each domain",sep=" ")  
     , xlim=c(0,50)
     ,las=1
     ,main= "")

f1=read.table("gencode.vM20.annotation_geneMerged.bed_count_in_T8_2Strand_p0.05_effect_imprinting.bed_cluster",header=F)
f1_0=read.table("T8_2Strand_p0.05_effect_imprinting.bed_cluster",header=F)

h=hist(c(rep(0,dim(f1_0)[1] - dim(f1)[1]),f1$V1),col="red" 
       ,density=25
       , breaks = seq(0,200,binSize)
       , freq = F
       ,ylab="Proportion of clusters"
       , xlim=c(0,50)
       ,xlab=paste("Number of","gencode gene annotations in each cluster",sep=" ")
       ,main= ""
       ,add=T
       ,las=1
       ,plot =FALSE
)
h$counts=h$counts/sum(h$counts)
plot(h,col="red" ,density=25, add=T)


legend("topright", 
       legend = c( "Imprinted","Strain effect"), 
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

#Organ_domain_counts
