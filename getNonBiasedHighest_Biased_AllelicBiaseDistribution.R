# function to merge multiple files based on common columns
multmerge = function(mypath, mypattern){
  filenames=list.files(path=mypath,  pattern = mypattern)
  datalist = lapply(filenames, function(x){read.table(x, head=T)})
  Reduce(function(x,y) {merge(x,y)}, datalist)}


## get binomial test p-value
load("../data-rpkms-001.RData")
pValue= multmerge(".", mypattern = glob2rx("*_pValue.txt"))
organs <- c("BN", "HT", "SK", "SP", "KD", "LV", "GI" ,"ST")

rpkm_s=rpkm[,grep("chr",colnames(rpkm))]
rpkm_s$BiasedOrgan <- rpkm$BiasedOrgan
rpkm_s$NonBH <- rpkm$NonBH
pValue=merge.data.frame(pValue, rpkm_s)

pValue$B.p=1 # the organ specific domain in Biased organ)
pValue$H.p=1 # the non-biased organ with highest expression level

for (i in 1:dim(pValue)[1]){
  pValue$B.p[i]=pValue[i,grep(pValue$BiasedOrgan[i],colnames(pValue))] #find the p value of biased organ
  pValue$H.p[i]=pValue[i,grep(pValue$NonBH[i],colnames(pValue))]#find the p value of non-biased organ with highest expression level
}

# make QQ plots
pdf("qqplot_pValue_Pooled_reads_1.pdf")
qqplot(-log(seq(0, 1, 1/5000),10), -log(pValue$H.p,10), col="red") 
abline(0,1)
dev.off()
pdf("hist_pValue_Pooled_reads_1.pdf")
hist(pValue$H.p, breaks = seq(0,1,1/100))
dev.off()


## get allele-specific reads
AS.Read= multmerge(".", mypattern = glob2rx("*_AlleleSpecificReads.txt"))
# get the Biased organ and the non-biased organ with highest expression level
AS.Read= merge(AS.Read, rpkm_s)

# identify the allelic bias states of the organ-specific blocks in the Biased organ (M or P) 
for (i in 1:dim(AS.Read)[1]){
  AS.Read$B.win[i] <- AS.Read[i,grep(paste(AS.Read$BiasedOrgan[i],"win",sep = "."),colnames(AS.Read))] 
}
AS.Read$B.win[AS.Read$B.win==1]="M"
AS.Read$B.win[AS.Read$B.win==2]="P"

# identify mat and pat reads of the non-biased highest organ and biased organ
for (i in 1:dim(AS.Read)[1]){
  AS.Read$B.mat[i] <-  AS.Read[i,grep(paste(AS.Read$BiasedOrgan[i],"mat",sep = "."),colnames(AS.Read))]
  AS.Read$B.pat[i] <-  AS.Read[i,grep(paste(AS.Read$BiasedOrgan[i],"pat",sep = "."),colnames(AS.Read))]
  AS.Read$H.mat[i] <-  AS.Read[i,grep(paste(AS.Read$NonBH[i],"mat",sep = "."),colnames(AS.Read))]
  AS.Read$H.pat[i] <-  AS.Read[i,grep(paste(AS.Read$NonBH[i],"pat",sep = "."),colnames(AS.Read))]
  }


# calulate effect size (mat/pat or pat/mat depends on m |p of the biased organ)
# mat/pat if the biased organ is M, otherwise pat/mat
for (i in 1:dim(AS.Read)[1]){
  if (AS.Read$B.win[i] == "M"){
    AS.Read$H.effectSize[i] = AS.Read$H.mat[i]/AS.Read$H.pat[i]
    AS.Read$B.effectSize[i] = AS.Read$B.mat[i]/AS.Read$B.pat[i]
    }
    else{
      AS.Read$H.effectSize[i] = AS.Read$H.pat[i]/AS.Read$H.mat[i]
      AS.Read$B.effectSize[i] = AS.Read$B.pat[i]/AS.Read$B.mat[i]
    }
}


AS.Read$H.matReads.TotalRatio = AS.Read$H.mat/(AS.Read$H.pat+AS.Read$H.mat)
AS.Read$B.matReads.TotalRatio = AS.Read$B.mat/(AS.Read$B.pat+AS.Read$B.mat)

#boxplot(log2(AS.Read$H.effectSize), log2(AS.Read$B.effectSize))
#abline(h=0)

# exmaine the effect size.
pdf("AllelicBias_effectSize_NonBiaseH_Biased.pdf")
par(mar=c(6.1, 7.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
par(mgp=c(3,1,0))
par(cex.lab=2.2, cex.axis=2.2)
hist(log2(AS.Read$H.effectSize),  col = "blue", breaks = seq(-10,20,1/10), xlim = c(-10,10), xlab="Effect Size (log2)", main="")
hist(log2(AS.Read$B.effectSize), xlim = c(-10,10), col="red", density=25, add=T, breaks = seq(-10,20,1/10))
abline(v=0, col="yellow", lty=3, lwd=2)
legend("topleft", 
       legend = c( "Biased","NonBiased"), 
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

pdf("AllelicBias_matRatio_NonBiaseH_Biased.pdf")
par(mar=c(6.1, 7.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
par(mgp=c(3,1,0))
par(cex.lab=2.2, cex.axis=2.2)
hist(AS.Read$H.matReads.TotalRatio, breaks = seq(0,1,1/30),col="blue", xlab="mat read ratio", main="")
hist(AS.Read$B.matReads.TotalRatio, breaks = seq(0,1,1/30), col = "red", density = 25, add=T)
legend("topright", 
       legend = c( "Biased","NonBiased"), 
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