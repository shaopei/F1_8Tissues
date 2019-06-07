args=(commandArgs(TRUE))
setwd(args[1])
df_fp=args[2] #"BN_MB6_paires_within_cluster1M_pairsWithin1M.txt"
maxD=as.numeric(args[3]) #3000
pdf_name=args[4] # TSS within maxD distance
pdf_name2=args[5] #run over



#setwd("~/Box Sync/Danko_lab_work/F1_8Tissues/Tunit_pair_analysis")
#df_fp="BN_MB6_paires_within_cluster1M_pairsWithin1M.txt"
print(df_fp)

df=read.table(df_fp, header = T)
head(df)

#check from same cluster
if (sum(df$ClusterID1!=df$ClusterID2) != 0) {   #should be 0
  warning("Paires Should from same cluster")
}


#check from different tunit
x=c(levels(df$Name1), levels(df$Name2))
x=x[!duplicated(x)]
df$Name1 <- factor(df$Name1, levels=x)
df$Name2 <- factor(df$Name2, levels=x)
if (sum(df$Name1 ==df$Name2) != 0) {   #should be 0
  warning("Paires Should from different tunit")
}


df=df[df$TSS_distance <= maxD,]
#View(subdf)

pdf(pdf_name)
d=30 #30bp
print(d)
print(maxD)
breaks = seq(0,maxD+d,d)
print (breaks)
hist(df$TSS_distance[df$CoDiS_group=="Con"], main=df_fp,xlab = "distance between TSS of Tunit pairs within 1M", 
     breaks = breaks , freq = T, col="blue", density=25,angle=135, add=F)
hist(df$TSS_distance[df$CoDiS_group=="Dis"],  main=df_fp,xlab = "distance between TSS of Tunit pairs within 1M", 
     breaks = breaks , freq = T, col="dark orange", density=25,angle=45, add=T)

#hist(df$TSS_distance[df$CoDiS_group=="OneS"],  breaks = breaks , freq = F, add=T)


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


df=df[df$TSS_runover=="True",]
pdf(pdf_name2)
hist(df$TSS_distance[df$CoDiS_group=="Con"], main=df_fp,xlab = "distance between TSS of runover Tunit pairs within 1M", 
     breaks = breaks , freq = T, col="blue", density=25,angle=135, add=F)
hist(df$TSS_distance[df$CoDiS_group=="Dis"],  main=df_fp,xlab = "distance between TSS of runover Tunit pairs within 1M", 
     breaks = breaks , freq = T, col="dark orange", density=25,angle=45, add=T)

#hist(df$TSS_distance[df$CoDiS_group=="OneS"],  breaks = breaks , freq = F, add=T)


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