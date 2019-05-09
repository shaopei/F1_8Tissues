setwd("~/Box Sync/Danko_lab_work/F1_8Tissues/UpSetR")
fp="T8_closest_io_minus_p0.05_effect_imprinting.bed"
df=read.table(fp)


hist(log10(df$V9[(df$V9>0)]))



Uplimit=5e06
hist(f1$V1, col="blue",
     breaks = seq(0,Uplimit,25000), 
     freq = FALSE,
     xlab="distance(bp)",
     main="distance to nearest region") 
Uplimit=1e05
Uplimit=1000000
hist(f1$V1[f1$V1<= Uplimit], col="blue",
     breaks = seq(0,Uplimit,5000), 
     freq = FALSE,
     xlab="distance(bp)",
     main="distance to nearest region")