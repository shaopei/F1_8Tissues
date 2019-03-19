# R --vanilla --slave --args ${tmp}_plus.bedGraph ${tmp}_minus.bedGraph readCount.${tmp} <ReadCounts_bedGraph.R
args=(commandArgs(TRUE))

a=read.table(args[1], sep="\t")
b=read.table(args[2], sep="\t")
cat(sum(a$V4)-sum(b$V4),file=args[3])

