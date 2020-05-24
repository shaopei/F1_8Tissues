import numpy as np
from math import *
from sys import argv
from multiprocessing import Pool

infp=argv[1]
outfp=argv[2]

data=np.array([l.strip().split("\t") for l in open(infp).readlines()])

chrom_list= list(set(data[:,0]))
chrom_list.sort()
data_dic={}
for c in chrom_list:
    data_dic[c] =  data[data[:,0]==c,:]
data=[]

def RemoveTSNwithMaxAllelicDifference(data_key):
    result=[]
    data=data_dic[data_key]
    count=1
    while data.shape[0] > 0:
        chrom, chromStart, chromEnd = data[0][0:3]
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
        
        diff= abs(subdata[:,7].astype(int) - subdata[:,15].astype(int))
        # remove the row with max allelic difference. If there are multiple base with the same diff, only first one will be removed
        new_subdata = np.delete(subdata, np.where(diff == np.amax(diff))[0][0], 0)
        for l in new_subdata:
            result.append("\t".join(l))
    return result

pool = Pool(processes=25)
pool_output = pool.map(RemoveTSNwithMaxAllelicDifference, chrom_list)
# pool_output looks like [[t, new_T_list,new_P_list, p_Y_f_list],...]
pool.close() # no more tasks
pool.join()



with open(outfp, "w") as out:
    for line in pool_output:
        out.write("\n".join(line))
        out.write("\n")
