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
     ,right = FALSE
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
     ,right = FALSE
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
       ,right = FALSE
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
       ,right = FALSE
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
#library(UpSetR)
setwd("~/Box Sync/Danko_lab_work/F1_8Tissues/UpSetR")
par(mar=c(6.1, 7.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
par(mgp=c(3,1,0))
par(cex.lab=2.2, cex.axis=2.2)
# imprinting cluster regardless of strandness
Tissue_list=c( "BN","SP","HT","SK","KD","ST","GI","LV")
df=read.table("T8_2Strand_p0.05_effect_imprinting.bed_cluster", header = F)

for (kkk in Tissue_list){
  df$tmp=0
  df$tmp[grepl(kkk, df$V6)]=1
  colnames(df)[grep("tmp", colnames(df))]=kkk
}
#upset(df, nsets = 8, sets =Tissue_list, keep.order = T, order.by = "degree")
barplot(colSums(df[ , match(Tissue_list , names(df) ) ] ), col="dark red", las=1, xlab= "Organ", ylab="Domain counts"  )

df$TissueCounts= rowSums(df[ , match(Tissue_list , names(df) ) ] )  
h1=hist(df$TissueCounts[df$TissueCounts >=1], breaks = seq(1,9,1),
 right=FALSE5,plot =FALSE)

h1$counts=h1$counts/sum(h1$counts)

h1_sub=hist(df$TissueCounts[df$TissueCounts >1], breaks = seq(2,9,1),
          right=FALSE,plot =FALSE)
h1_sub$counts=h1_sub$counts/sum(h1_sub$counts)

# strain effect cluster regardless of strandness
Tissue_list=c( "BN","SP","HT","SK","KD","ST","GI","LV")
df=read.table("T8_2Strand_p0.05_effect_strain.bed_cluster", header = F)
for (kkk in Tissue_list){
  df$tmp=0
  df$tmp[grepl(kkk, df$V6)]=1
  colnames(df)[grep("tmp", colnames(df))]=kkk
}
#upset(df, nsets = 8, sets =Tissue_list, keep.order = T, order.by = "degree", nintersects=100)
barplot(colSums(df[ , match(Tissue_list , names(df) ) ] ), col = "blue", las=1, xlab= "Organ", ylab="Domain counts"  )
df$TissueCounts= rowSums(df[ , match(Tissue_list , names(df) ) ] )  

h2=hist(df$TissueCounts[df$TissueCounts >=1], breaks = seq(1,9,1),
        right=FALSE,plot =FALSE)
h2$counts=h2$counts/sum(h2$counts)
h2_sub=hist(df$TissueCounts[df$TissueCounts >1], breaks = seq(2,9,1),
            right=FALSE,plot =FALSE)
h2_sub$counts=h2_sub$counts/sum(h2_sub$counts)

## Strain_imprinted_domains
# Number of organs with allelic biased blocks in the domain
plot(h2,col="blue", xlab="Number of organs with allelic biased blocks in the domain", ylab="Proportion of clusters", main="",
     las=1)
plot(h1,col="red", density=25, add=T)

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
# Number of organs with allelic biased blocks in the domain (at least two organs)
plot(h2_sub,col="blue", xlab="Number of organs with allelic biased blocks in the domain", ylab="Proportion of clusters", main="",
     las=1)
plot(h1_sub,col="red", density=25, add=T)

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
