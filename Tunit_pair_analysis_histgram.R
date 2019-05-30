args=(commandArgs(TRUE))
setwd(args[1])
df_fp=args[2] #"BN_MB6_paires_within_cluster500K.txt"
pdf_name=args[3]

#df=read.table("BN_MB6_paires_within_cluster500K.txt")
df=read.table(df_fp)

print(df_fp)
print(head(df))
# V1       V2       V3                V4 V5 V6  V7   V8       V9      V10
# 1 chr3 66968150 66981450 preds_minus_14591  S  - 344 chr3 66981350 67371050
# 2 chr3 66968150 66981450 preds_minus_14591  S  - 344 chr3 67373900 67374000
# 3 chr3 66968150 66981450 preds_minus_14591  S  - 344 chr3 67373950 67402500
# 4 chr3 66968150 66981450 preds_minus_14591  S  - 344 chr3 67427050 67430050
# 5 chr3 66968150 66981450 preds_minus_14591  S  - 344 chr3 67430000 67430200
# 6 chr3 66968150 66981450 preds_minus_14591  S  - 344 chr3 67430250 67464050
# 
# V11 V12 V13 V14  V15    V16
# 1  preds_plus_14593   S   + 344 OneS    100
# 2 preds_minus_14592   S   - 344 OneS 392550
# 3  preds_plus_14594   S   + 344 OneS 392500
# 4 preds_minus_14593   S   - 344 OneS 448600
# 5  preds_plus_14595   S   + 344 OneS 448550
# 6  preds_plus_14596   S   + 344 OneS 448800





#check from same cluster
if (sum(df$V7!=df$V14) != 0) {   #should be 0
  warning("Paires Should from same cluster")
}


#check from different tunit
x=c(levels(df$V11), levels(df$V4))
x=x[!duplicated(x)]
df$V11 <- factor(df$V11, levels=x)
df$V4 <- factor(df$V4, levels=x)
if (sum(df$V4 ==df$V11) != 0) {   #should be 0
  warning("Paires Should from different tunit")
}

#histogram of the distance df$V16

pdf(pdf_name)
d=1000000 #1M
breaks = seq(0,max(df$V16)+d,d)
hist(df$V16[df$V15=="Con"], main=df_fp,xlab = "distance between TSS of Tunit pairs within cluster", 
     breaks = breaks , freq = F, col="blue", density=25,angle=135)
hist(df$V16[df$V15=="Dis"],  breaks = breaks , freq = F, col="dark orange", density=25,angle=45, add=T)
hist(df$V16[df$V15=="OneS"],  breaks = breaks , freq = F, add=T)


legend("topright", 
       legend = c("OneS", "Concordant", "Discordant"),
       #pch=c(15,15),
       cex=2, 
       lty=c(0,0),
       #bty="n",
       lwd=1.5, 
       density=c(25,25,25),
       angle=c(180, 135,45),
       #angle=45,
       fill=c("white","blue", "dark orange")
       , bty = "n"
)
dev.off()

t=df


maxD=3000000
df=t[t$V16<=maxD,]
pdf(paste(pdf_name, "_PairesWithin3Mb.pdf", sep=""))
d=100000 #0.1M
#d=50000
breaks = seq(0,maxD,d)
hist(df$V16[df$V15=="Con"], main=df_fp,xlab = "distance between TSS of Tunit pairs within cluster", 
     breaks = breaks ,
     freq = F, col="blue", density=25,angle=135)
hist(df$V16[df$V15=="Dis"] , breaks = breaks ,freq = F, col="dark orange", density=25,angle=45, add=T)
hist(df$V16[df$V15=="OneS"], breaks = breaks , freq = F, add=T)

 
legend("topright", 
       legend = c("OneS", "Concordant", "Discordant"),
       #pch=c(15,15),
       cex=2, 
       lty=c(0,0),
       #bty="n",
       lwd=1.5, 
       density=c(25,25,25),
       angle=c(180, 135,45),
       #angle=45,
       fill=c("white","blue", "dark orange")
       , bty = "n"
)
dev.off()


maxD=1000000
df=t[t$V16<=maxD,]
pdf(paste(pdf_name, "_PairesWithin1Mb.pdf", sep=""))
d=50000 #50K
breaks = seq(0,maxD,d)
hist(df$V16[df$V15=="Con"], main=df_fp,xlab = "distance between TSS of Tunit pairs within cluster", 
     breaks = breaks ,
     freq = F, col="blue", density=25,angle=135)
hist(df$V16[df$V15=="Dis"] , breaks = breaks ,freq = F, col="dark orange", density=25,angle=45, add=T)
hist(df$V16[df$V15=="OneS"], breaks = breaks , freq = F, add=T)


legend("topright", 
       legend = c("OneS", "Concordant", "Discordant"),
       #pch=c(15,15),
       cex=2, 
       lty=c(0,0),
       #bty="n",
       lwd=1.5, 
       density=c(25,25,25),
       angle=c(180, 135,45),
       #angle=45,
       fill=c("white","blue", "dark orange")
       , bty = "n"
)
dev.off()



