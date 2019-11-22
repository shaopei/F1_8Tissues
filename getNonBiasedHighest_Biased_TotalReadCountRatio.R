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

hist(log2(rpkm$NonBH_B_ratio), xlab="log2(rpkm non-biased organ with highest expression level/Biased organ)", breaks = seq(-10,10,1/10))
dev.off()