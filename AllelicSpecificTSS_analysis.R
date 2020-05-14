setwd("~/Box Sync/Danko_lab_work/F1_8Tissues/Initiation/")
tss=read.table(file = "BN_allReads_TSS_5mat5pat_uniq_WithAsMaxTSNunionAsTSSfdr0.1_maskedVSunmasked_pValue.bed")
colnames(tss)[7]="masked_p_value"
colnames(tss)[9]="unmasked_p_value"
colnames(tss)[8]="masked_fdr"
colnames(tss)[10]="unmasked_fdr"
plot( -log10(tss$unmasked_p_value),-log10(tss$masked_p_value), pch=19, col=rgb(0,0,0,alpha = 0.25))
abline(v=-log10(max(tss$unmasked_p_value[tss$unmasked_fdr<=0.1])), col="red")
abline(0,1, col="blue")

p_value_cut = max(tss$unmasked_p_value[tss$unmasked_fdr<=0.1])  


plot( (tss$unmasked_p_value),(tss$masked_p_value), pch=19, col=rgb(0,0,0,alpha = 0.25))
abline(v=max(tss$unmasked_p_value[tss$unmasked_fdr<=0.1]), col="red")
abline(0,1, col="blue")

plot( -log10(tss$unmasked_p_value[tss$unmasked_fdr<=0.1]),-log10(tss$masked_p_value[tss$unmasked_fdr<=0.1]), pch=19, col=rgb(0,0,0,alpha = 0.25))
abline(0,1, col="blue")

# how many As.TSS became not significant using KS test after MASK one bp with largest difference? 685
# Not AS.TSS after masked 1 bp, driven by a single base
hist(tss$unmasked_p_value, col="blue", breaks = seq(0,1,0.01))
hist(tss$masked_p_value, add=T, col="red", density = 25, breaks = seq(0,1,0.01))
legend("topright", 
       legend = c("unmaked", "masked"),
       #pch=c(15,15),
       cex=2, 
       lty=c(0,0),
       #bty="n",
       lwd=1.5, 
       density=c(10000,25),
       angle=c(180,45),
       #angle=45,
       fill=c("blue","red")
       , bty = "n"
)

sub_tss=tss[(tss$masked_p_value>0.01 & tss$unmasked_fdr<=0.1),]
dim(sub_tss)
sub_tss=sub_tss[order(sub_tss$masked_p_value, decreasing = TRUE),]

# STILL AS.TSS after masked 1 bp, driven by more than one bases

still_tss=tss[(tss$masked_p_value<=p_value_cut),]
dim(still_tss)
dim(tss[(tss$masked_p_value<=p_value_cut & tss$unmasked_fdr <=0.1),])
dim(tss[(tss$masked_p_value>p_value_cut & tss$unmasked_fdr <=0.1),])

# how many noAs.TSS.with_as.maxTSN became as.TSS after mask? 18 out of 324
# Gain AS.TSS after mask one bp with biggest allelic difference
hist(tss$unmasked_p_value[tss$unmasked_fdr>0.1], col="blue", breaks = seq(0,1,0.01))
hist(tss$masked_p_value[tss$unmasked_fdr>0.1], add=T, col="red", density = 25, breaks = seq(0,1,0.01))

dim(tss[(tss$unmasked_fdr>0.1),])
dim(tss[(tss$masked_p_value<=p_value_cut & tss$unmasked_fdr>0.1),])
dim(tss[(tss$masked_p_value >p_value_cut & tss$unmasked_fdr>0.1),])
dim(tss[(tss$unmasked_fdr<=0.1),])


f = tss[(tss$masked_p_value > p_value_cut),]
f=f[order(f$masked_p_value, decreasing = FALSE),]
View(f)
dim(tss[(tss$masked_p_value >  0.05),])

hist(tss$masked_p_value, col="red", density = 25, breaks = seq(0,1,0.01), ylim=c(0,200))
