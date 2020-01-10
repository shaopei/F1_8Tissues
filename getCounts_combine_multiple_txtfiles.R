multmerge = function(mypath, mypattern){
  filenames=list.files(path=mypath,  pattern = mypattern)
  datalist = lapply(filenames, function(x){read.table(x, head=T)})
  Reduce(function(x,y) {merge(x,y)}, datalist)}

df= multmerge(".", mypattern = glob2rx("*_readcounts.txt"))


organs <- c("BN", "HT", "SK", "SP", "KD", "LV", "GI" ,"ST")

## Gets counts
# combine counts from plus amd minus strand
# combine counts from PB6 and MB6
counts <- df[,1:4]
for(o in organs) {
  counts$temp = rowSums(df[,grep(o, colnames(df))])
  colnames(counts)[grep("temp", colnames(counts))] = o
 # counts <- cbind(counts, countBigWig(f, bodies, rpkm=FALSE))
}
#save.image("data-counts.RData")
#load("data-counts.RData")

## Gets RPKMs
rpkm  <-counts
targetCol = grep( "chr",colnames(counts), invert = T)
rpkm[,targetCol] <- counts[,targetCol] * (1000/(counts$chrmEnd - counts$chrmStart)) * (1e6/colSums(counts[targetCol]))

osdBlock = read.table("T8_2Strand_p0.05_effect_strain_withStrandness.bed_organSpecific", header = F)
colnames(osdBlock)=c("chrm", "chrmStart", "chrmEnd", "BiasedOrgan", "_", "chrStrand")
rpkm <- merge(rpkm, osdBlock[,grep( "_",colnames(osdBlock), invert = T)])

rpkm$B=0 # rpkm of the organ specific domain in Biased organ)
rpkm$H=0 # rpkm of the non-biased organ with highest expression level
rpkm$NonBH = "A"
for (i in 1:dim(rpkm)[1]){
  rpkm$B[i]=rpkm[i,grep(rpkm$BiasedOrgan[i],colnames(rpkm))]
  rpkm$H[i]=max(rpkm[,targetCol][i,grep(rpkm$BiasedOrgan[i],colnames(rpkm)[targetCol], invert = T )])
  rpkm$NonBH[i]=colnames(rpkm)[targetCol][grep(rpkm$H[i],rpkm[i,targetCol])]
}

save.image("data-rpkms-001.RData")
#load("data-rpkms.RData")
rpkm$NonBH_B_ratio = rpkm$H/rpkm$B
pdf("allellicBias_expression_level3.pdf")
par(mar=c(6.1, 7.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
par(mgp=c(3,1,0))
par(cex.lab=2.2, cex.axis=2.2)

hist(log2(rpkm$NonBH_B_ratio)#, xlab="log2(rpkm non-biased organ with highest expression level/Biased organ)"
     , breaks = seq(-10,10,1/10), xlab="log2(rpkm NonBiasHighest/Biased)", main="")
abline(v=0, col="red")
dev.off()

## get binomial test p-value
pValue= multmerge(".", mypattern = glob2rx("*_pValue.txt"))
organs <- c("BN", "HT", "SK", "SP", "KD", "LV", "GI" ,"ST")

# use min p-value from plus and minus strand, from PB6 and MB6
#pValue <- df[,1:3]
#for(o in organs) {
#  pValue$temp = apply(df[,grep(o, colnames(df))], 1, FUN=min) # use min p-value from plus and minus strand, from PB6 and MB6
#  colnames(pValue)[grep("temp", colnames(pValue))] = o
# }
#save.image("data-pValue.RData")
load("data-rpkms-001.RData")
load("data-pValue.RData")
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


## 
AS.Read= multmerge(".", mypattern = glob2rx("*_AlleleSpecificReads.txt"))
save.image("data-AS.Read.RData")
load("data-AS.Read.RData")
AS.Read= merge(AS.Read, rpkm_s)

# identify M or P of the organ-specific blocks
for (i in 1:dim(AS.Read)[1]){
  AS.Read$B.win[i] <- AS.Read[i,grep(paste(AS.Read$BiasedOrgan[i],"win",sep = "."),colnames(AS.Read))] #find the p value of biased organ
  #pValue$H.p[i]=pValue[i,grep(pValue$NonBH[i],colnames(pValue))]#find the p value of non-biased organ with highest expression level
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