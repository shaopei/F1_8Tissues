import numpy as np
from sys import argv

def bed_overlap(bed1, bed2):
    s1, e1 = bed1
    s2, e2 = bed2
    s1 = int(s1)
    s2 = int(s2)
    e1 = int(e1)
    e2 = int(e2)
    if (s2 < e1 and e1 <= e2) or (s2 <= s1 and s1 < e2):
        return True
    else: 
        s1, e1 = bed2
        s2, e2 = bed1
        s1 = int(s1)
        s2 = int(s2)
        e1 = int(e1)
        e2 = int(e2)
        return (s2 < e1 and e1 <= e2) or (s2 <= s1 and s1 < e2)



def bed_distance(bed1, bed2):
    s1, e1 = bed1
    s2, e2 = bed2
    if bed_overlap(bed1, bed2):
        return 0
    else:
        return min(abs(int(s1)-int(e2)), abs(int(s2)-int(e1))) +1

def TSS_distance(Sbed1, Sbed2):
    s1, e1, strand1 = Sbed1
    s2, e2, strand2 = Sbed2
    s1 = int(s1)
    s2 = int(s2)
    e1 = int(e1)
    e2 = int(e2)
    if strand1 == "+" and strand2 == "+":
        return abs(s1-s2)
    elif strand1 == "-" and strand2 == "-":
        return abs(e1-e2)
    elif strand1 == "+" and strand2 == "-":
        return abs(s1-e2)
    elif strand1 == "-" and strand2 == "+":
        return abs(e1-s2)
    else:
        return "ERROR in TSS_distance"


def Sbed1_TSS_run_over_by_Sbed2(Sbed1, Sbed2): #Strand Bed1, Strand Bed2
    s1, e1, strand1 = Sbed1
    s2, e2, strand2 = Sbed2
    s1 = int(s1)
    s2 = int(s2)
    e1 = int(e1)
    e2 = int(e2)
    if bed_overlap((s1,e1), (s2,e2)):
        if strand1 == "+" and strand2 == "-":
            return s2 <= s1
        elif strand1 == "-" and strand2 == "+":
            return e1<= e2
        else:
            return False
    else:
        return False



infp=argv[1] #"BN_PB6_paires_within_cluster1M.txt"
outfp=argv[2] #"BN_PB6_paires_within_cluster1M_overlapped_pairs.txt"
 

data=np.array([l.strip().split("\t") for l in open(infp).readlines()])
data.shape
data=data[data[:,14]!="TwoS"] # remove lines with TwoS
data.shape
data=data[data[:,5]!=data[:,12]] # keep line with tunit on different strand (+,-)
data.shape
data=data[data[:,0]!="chrX"] # remove lines on chrX
data.shape

"""
#pairs overlap
TSS_d=[]
TSS_run_over=[]
CDS_group=[]
newData=[]
bed_d=[]

for l in data:
    if bed_overlap(l[1:3], l[8:10]):
        if l[4]=="S": # for OneS, Keep S as Sbed1 
            Sbed1=(l[8], l[9], l[12])
            Sbed2=(l[1], l[2], l[5])
        else:
            Sbed1=(l[1], l[2], l[5])
            Sbed2=(l[8], l[9], l[12])
        newData.append(l)
        CDS_group.append(l[14])
        TSS_d.append(TSS_distance(Sbed1, Sbed2))
        bed_d.append(bed_distance(Sbed1[0:2], Sbed2[0:2]))
        TSS_run_over.append(Sbed1_TSS_run_over_by_Sbed2(Sbed1, Sbed2))



td = np.array(TSS_d)
tr = np.array(TSS_run_over)
cds  = np.array(CDS_group)
for g in ["Con","Dis", "OneS"]:
    print g, sum(data[:,14]==g), "Overlap", sum(cds==g), round(sum(cds==g)*1.0/sum(data[:,14]==g),4), "RunOver", sum(tr[cds==g]), round(sum(tr[cds==g])*1.0/(sum(cds==g)),2)

with open(outfp, 'w') as out:
    out.write('\t'.join(["chrm1", "chrmStart1", "chrmEnd1", "Name1", "WinP1","Strand1", "ClusterID1","chrm2", "chrmStart2", "chrmEnd2", "Name2", "WinP2","Strand2", "ClusterID2", "CoDiS_group", "TSS_distance", "bed_distance", "TSS_runover"]) )
    out.write('\n')
    for i in xrange(len(newData)):
        out.write('\t'.join(list(newData[i])+[str(bed_d[i]),str(TSS_run_over[i])])) 
        out.write('\n')
"""

# Pairs Within 1Mb
TSS_d=[]
TSS_run_over=[]
CDS_group=[]
newData=[]
bed_d=[]

for l in data:
    if bed_distance(l[1:3], l[8:10]) <= 1000000:
        if l[4]=="S": # for OneS, Keep S as Sbed2
            Sbed1=(l[8], l[9], l[12])
            Sbed2=(l[1], l[2], l[5])
        else:
            Sbed1=(l[1], l[2], l[5])
            Sbed2=(l[8], l[9], l[12])
        #print l
        #print Sbed1, Sbed2
        #print "TSS distance", TSS_distance(Sbed1, Sbed2)
        #print "bed distance", bed_distance(Sbed1[0:2], Sbed2[0:2])
        newData.append(l)
        CDS_group.append(l[14])
        TSS_d.append(TSS_distance(Sbed1, Sbed2))
        bed_d.append(bed_distance(Sbed1[0:2], Sbed2[0:2]))
        TSS_run_over.append(Sbed1_TSS_run_over_by_Sbed2(Sbed1, Sbed2))

import scipy.stats as stats


td = np.array(TSS_d)
tr = np.array(TSS_run_over)
cds  = np.array(CDS_group)
for g in ["Con","Dis", "OneS"]:
    if sum(data[:,14]==g)>0 :
        print g, sum(data[:,14]==g), "WithIn1M", sum(cds==g), "RunOver", sum(tr[cds==g])

oddsratio, pvalue = stats.fisher_exact([[sum(cds=="Con"), sum(tr[cds=="Con"])], [sum(cds=="Dis"), sum(tr[cds=="Dis"])]])
print oddsratio, pvalue



with open(outfp, 'w') as out:
    out.write('\t'.join(["chrm1", "chrmStart1", "chrmEnd1", "Name1", "WinP1","Strand1", "ClusterID1","chrm2", "chrmStart2", "chrmEnd2", "Name2", "WinP2","Strand2", "ClusterID2", "CoDiS_group", "TSS_distance", "bed_distance", "TSS_runover"]) )
    out.write('\n')
    for i in xrange(len(newData)):
        out.write('\t'.join(list(newData[i])+[str(bed_d[i]),str(TSS_run_over[i])])) 
        out.write('\n')


