# output the baseS (can be more than 1) with max reads count within each TSS (define by the infout bed6 file chrom, chromStart, chromEnd, _, _, strand)

import numpy as np
from math import *
from sys import argv
infp=argv[1]
outfp=argv[2]


data=np.array([l.strip().split("\t") for l in open(infp).readlines()])
tunit_list= list(set(data[:,9]))
data_dic={}
for t in tunit_list:
    data_dic[t] =  data[data[:,9]==t,:]

data=[]
    

with open(outfp, "w") as out:
    for t in tunit_list:
        #result=[]
        data=data_dic[t]
        #count=1
        readCount = data[:,3].astype(int)
        maxTSNs=data[readCount==max(readCount)]
        if len(maxTSNs) > 0 :
            for s in maxTSNs:
                out.write("\t".join(s))
                out.write("\n")
