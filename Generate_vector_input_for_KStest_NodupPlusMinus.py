import numpy as np
from math import *
from sys import argv
from multiprocessing import Pool
infp=argv[1]
outfp=argv[2]
span=100

sample="_".join(infp.split("_")[0:2])
data=np.array([l.strip().split("\t") for l in open(infp).readlines()])
chrom_list= list(set(data[:,0]))
chrom_list.sort()
data_dic={}
for c in chrom_list:
    data_dic[c] =  data[data[:,0]==c,:]

data=[]

def getVector(data_key):
    result=[]
    data=data_dic[data_key]
    count=1
    while data.shape[0] > 0:
        chrom, chromStart, chromEnd, name = data[0][0:4]
        #c0 = data[:,0]==chrom
        c1 = data[:,1]==chromStart
        c2 = data[:,2]==chromEnd
        #name = chrom+"_"+str(count)
        inside = np.logical_and(c1,c2)
        outside = np.invert(inside)
        subdata = data[inside,:]
        data = data[outside,:]
        sub_plus = subdata[subdata[:,5] == "+", :]
        sub_minus = subdata[subdata[:,5] == "-", :]
        
        sp=[]
        for s in sub_plus:
            sp = sp + [s[6]]*int(s[7])
        
        sm=[]
        for s in sub_minus:
            sm = sm + [s[6]]*int(s[7])
        
        if len(sm) > 0:
            result.append("\t".join([chrom, chromStart, chromEnd, name,".", "-", ",".join(sm)]))
        else:
            result.append("\t".join([chrom, chromStart, chromEnd, name,".", "-", "0"]))
        if len(sp) > 0:
            result.append("\t".join([chrom, chromStart, chromEnd, name,".", "+", ",".join(sp)]))
        else:
            result.append("\t".join([chrom, chromStart, chromEnd, name,".", "+", "0"]))
        count = count + 1
    return result

def getVector_noPM_duplicate(data_key):
    result=[]
    data=data_dic[data_key]
    count=1
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
        
        sp=[]
        for s in subdata:
            sp = sp + [s[6]]*int(s[7])
        
        if len(sp) > 0:
            result.append("\t".join([chrom, chromStart, chromEnd, name,".", strand, ",".join(sp)]))
        else:
            result.append("\t".join([chrom, chromStart, chromEnd, name,".", strand, "0"]))
        count = count + 1
    return result

pool = Pool(processes=25)
pool_output = pool.map(getVector_noPM_duplicate, chrom_list)
# pool_output looks like [[t, new_T_list,new_P_list, p_Y_f_list],...]
pool.close() # no more tasks
pool.join()



with open(outfp, "w") as out:
    for line in pool_output:
        out.write("\n".join(line))
        out.write("\n")
