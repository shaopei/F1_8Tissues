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

def getSpan(data_key):
    result=[]
    data=data_dic[data_key]
    count=1
    while data.shape[0] > 0:
        chrom, chromStart, chromEnd = data[0][0:3]
        #c0 = data[:,0]==chrom
        c1 = data[:,1]==chromStart
        c2 = data[:,2]==chromEnd
        name = chrom+"_"+str(count)
        inside = np.logical_and(c1,c2)
        outside = np.invert(inside)
        subdata = data[inside,:]
        data = data[outside,:]
        #sub_plus = subdata[subdata[:,5] == "+", :]
        #sub_minus = subdata[subdata[:,5] == "-", :]
        #m = (int(sub_plus[:,-2][np.argmax(sub_plus[:,-1].astype(int))]) + int(sub_minus[:,-2][np.argmax(sub_minus[:,-1].astype(int))]))/2
        #result.append("\t".join([chrom, str(int(chromStart)+m-span), str(int(chromStart)+m), name,".", "-"]))
        #result.append("\t".join([chrom, str(int(chromStart)+m), str(int(chromStart)+m+span), name, ".", "+"]))
        result.append("\t".join([chrom, chromStart, chromEnd, name,".", "-"]))
        result.append("\t".join([chrom, chromStart, chromEnd, name,".", "+"]))
        count = count + 1
    return result


pool = Pool(processes=50)
pool_output = pool.map(getSpan, chrom_list)
# pool_output looks like [[t, new_T_list,new_P_list, p_Y_f_list],...]
pool.close() # no more tasks
pool.join()



with open(outfp, "w") as out:
    for line in pool_output:
        out.write("\n".join(line))
        out.write("\n")
