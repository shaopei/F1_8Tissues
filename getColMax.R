#R --vanilla --slave --args file col < getColMax.R &
args=(commandArgs(TRUE))
f=args[1]
c=as.integer(args[2])


df=read.table(f, sep="\t")
cat(max(df[,c]))
