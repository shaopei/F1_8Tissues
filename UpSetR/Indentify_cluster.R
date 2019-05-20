setwd("~/Box Sync/Danko_lab_work/F1_8Tissues/UpSetR")
fp="T8_closest_io_minus_p0.05_effect_imprinting.bed"
mdf=read.table(fp)
hist(log10(mdf$V9[(mdf$V9>0)]))

fp="T8_closest_io_plus_p0.05_effect_imprinting.bed"
pdf=read.table(fp)
hist(log10(pdf$V9[(pdf$V9>0)]))

d=c(pdf$V9, mdf$V9)
hist(log10(d[(d>0)]))
