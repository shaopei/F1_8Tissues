#
library(bigWig)

args <- commandArgs(trailingOnly=TRUE)

bw = load.bigWig(args[1])

for (chrom in bw$chroms)
    cat(chrom, "\n", sep='')

