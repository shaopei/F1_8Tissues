import numpy as np
from math import *
from sys import argv
from multiprocessing import Pool
infp=argv[1]
outfp=argv[2]


#sample="_".join(infp.split("_")[0:2])
data=np.array([l.strip().split("\t") for l in open(infp).readlines()])
chrom_list= list(set(data[:,0]))
chrom_list.sort()
data_dic={}
for c in chrom_list:
    data_dic[c] =  data[data[:,0]==c,:]

data=[]

def getVector_readLength(data_key):
    result=[]
    data=data_dic[data_key]
    count=0
    while data.shape[0] > 0:
        chrom, chromStart, chromEnd, name = data[0][0:4]
        strand=data[0][5]
        #c0 = data[:,0]==chrom
        c1 = data[:,1]==chromStart
        c2 = data[:,2]==chromEnd
        c5 = data[:,5]==strand
        #name = chrom+"_"+str(count)
        inside = np.logical_and(np.logical_and(c1,c2), c5)
        outside = np.invert(inside)
        subdata = data[inside,:]
        data = data[outside,:]
        
        read_length=[]
        for r in subdata[:,9]:
            read_length.append(str(len(r.split(",")[1])))
        
        read_length.sort()
        result.append("\t".join([chrom, chromStart, chromEnd, name,".", strand, ",".join(read_length)]))
        count = count + 1
    return result

pool = Pool(processes=25)
pool_output = pool.map(getVector_readLength, chrom_list)
# pool_output looks like [[t, new_T_list,new_P_list, p_Y_f_list],...]
pool.close() # no more tasks
pool.join()



with open(outfp, "w") as out:
    for line in pool_output:
        out.write("\n".join(line))
        out.write("\n")
