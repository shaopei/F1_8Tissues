#imput file format:
#chrm, chrmStart, chrmEnd, tunitName, Allele-bias, strand, clusterName.
#Seperate by tab. examlpes:
#chr1    3513700 3671800 preds_minus_1   S       -       1
#chr1    3671950 4459850 preds_plus_1    M       +       1
#chr1    3821300 4471500 preds_minus_2   S       -       1
#chr1    4487200 4493750 preds_minus_3   S       -       1
#chr1    4493750 4493850 preds_plus_2    S       +       1

# output file:
#chrm, chrmStart, chrmEnd, tunitName1, Allele-bias, strand, clusterName, chrm, chrmStart, chrmEnd, tunitName2, Allele-bias, strand, clusterName, Con/Dis/OneS, TSS_Distance

# python Pair_of_bedregions_withincluster.py input output

import numpy as np
from sys import argv

infp=argv[1] #"BN_all_h5.preds_2strands_MB6_ABconsistent_FisherMethodP0.05_cluster500K.bed"
outfp=argv[2] #"BN_paires_within_cluster.txt"
data=np.array([l.strip().split("\t") for l in open(infp).readlines()])
Column_Number=data.shape[1]
clustrer_site=7-1
tunit_name_site=4-1
tunit_AS_site=5-1


chrom = data[:,0]
bed_regions=data[:,1:3].astype(int)
cluster=data[:,clustrer_site]


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
        return (s2 < e1 and e1 <= e2) or (s2 <= s1 and s1 < e2)



def bed_distance(bed1, bed2):
    s1, e1 = bed1
    s2, e2 = bed2
    if bed_overlap(bed1, bed2):
        return 0
    else:
        return min(abs(int(s1)-int(e2)), abs(int(s2)-int(e1)))


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
    



# pair within each clusters
pair_within_cluster=[]
for c in set(cluster):
    data_c=data[cluster==c]
    bed_c = bed_regions[cluster == c]
    for i in xrange(len(data_c)):
        for j in xrange(i+1,len(data_c)):
            #print i, j
            pair_within_cluster.append(list(data_c[i])+list(data_c[j]))

p_list=pair_within_cluster
# Seperate into three group, Concordant (M,M)(P,P), Discordant (M,P)(P,M), OneS (S,M),(S,P)(P,S)(M,S)
pair_within_cluster=np.array(pair_within_cluster)
AS_group =  np.chararray(pair_within_cluster.shape[0],itemsize=4)
AS_group[:] = 'TwoS'

OneS = np.logical_or(pair_within_cluster[:,tunit_AS_site]!='S',pair_within_cluster[:,tunit_AS_site+Column_Number]!='S')
NoS = np.logical_and(pair_within_cluster[:,tunit_AS_site]!='S',pair_within_cluster[:,tunit_AS_site+Column_Number]!='S')
Con = np.logical_and(NoS, pair_within_cluster[:,tunit_AS_site]== pair_within_cluster[:,tunit_AS_site+Column_Number])
Dis = np.logical_and(NoS,pair_within_cluster[:,tunit_AS_site]!= pair_within_cluster[:,tunit_AS_site+Column_Number])
print "NoS", np.count_nonzero(NoS)
print "Con", np.count_nonzero(Con)
print "Dis", np.count_nonzero(Dis)

AS_group[OneS]="OneS"
AS_group[Con]="Con"
AS_group[Dis]="Dis"
print "OneS", np.count_nonzero(AS_group=="OneS")


# export the file
with open(outfp, 'w') as out:
    for i in xrange(len(p_list)):
        out.write('\t'.join(p_list[i]+[AS_group[i]]+[str(TSS_distance( p_list[i][1:3]+[p_list[i][5]],p_list[i][1+Column_Number:3+Column_Number]+[p_list[i][5+Column_Number]]))]))
#        out.write('\t'.join(p_list[i]+[AS_group[i]]+[str(bed_distance(p_list[i][1:3],p_list[i][1+Column_Number:3+Column_Number]))]))
        out.write('\n')

