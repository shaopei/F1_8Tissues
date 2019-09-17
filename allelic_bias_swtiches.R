setwd("~/Box Sync/Danko_lab_work/F1_8Tissues/")
df=read.table("allelic_bias_swtiches.txt", header = T)
for (i in 1:nrow(df)){
  a = fisher.test(matrix(c(df$im1[i],df$im2[i], df$st1[i], df$st2[i]), 2,2))
  df$p.value[i] = a$p.value
  df$odds.ratio=a$estimate
}