# output the baseS (can be more than 1) with max reads count within each TSS (define by the infout bed6 file chrom, chromStart, chromEnd, _, _, strand)

import numpy as np
from math import *
from sys import argv
infp=argv[1]
outfp=argv[2]


data=np.array([l.strip().split("\t") for l in open(infp).readlines()])
chrom_list= list(set(data[:,0]))
chrom_list.sort()
data_dic={}
for c in chrom_list:
    data_dic[c] =  data[data[:,0]==c,:]

data=[]
    

with open(outfp, "w") as out:
    for c in chrom_list:
        #result=[]
        data=data_dic[c]
        #count=1
        while data.shape[0] > 0:
            chrom, chromStart, chromEnd, _, _, strand = data[0][0:6]
            #strand=data[0][5]
            #c0 = data[:,0]==chrom
            c1 = data[:,1]==chromStart
            c2 = data[:,2]==chromEnd
            c5 = data[:,5]==strand
            #name = chrom+"_"+str(count)
            inside = np.logical_and(np.logical_and(c1, c2), c5)
            outside = np.invert(inside)
            subdata = data[inside,:]
            subdata_np = np.array(subdata)
            aveScore = subdata_np[:,18].astype(float)
            data = data[outside,:]
            #maxReadCounts= max(readCount)
            maxScore=subdata[aveScore==max(aveScore)]
            if len(maxScore) > 0 :
            #    for s in maxScore[0]:  # only use the first, if tie
                out.write("\t".join(maxScore[0]))
                out.write("\n")
