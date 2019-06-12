import numpy as np
from sys import argv
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

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


infp=argv[1] #"BN_MB6_paires_within_cluster1M.txt"
pdfName=argv[2]
#outfp=argv[2] #"BN_MB6_paires_within_cluster1M_overlapped_pairs.txt"
 
#reduce memory use
data=[]
with open(infp) as infile:
    for line in infile:
        l = line.strip().split("\t")
        #if l[14]!="TwoS": # remove lines with TwoS
        if l[5]!= l[12]: # keep line with tunit on different strand (+,-)
            if l[0]!="chrX": # remove lines on chrX
                data.append(l)


data=np.array(data)



#data=np.array([l.strip().split("\t") for l in open(infp).readlines()])
#data.shape
#data=data[data[:,14]!="TwoS"] # remove lines with TwoS
#data.shape
#data=data[data[:,5]!=data[:,12]] # keep line with tunit on different strand (+,-)
#data.shape
#data=data[data[:,0]!="chrX"] # remove lines on chrX
#data.shape
# ['chr7', '40773500', '40773900', 'preds_minus_20350', 'S', '-', '216',
#  'chr7', '41499100', '41525200', 'preds_plus_20310', 'P', '+', '216',
# 'OneS', '725200']


# make a copy for each pair in the reverser order. i.e. TunitA TunitB to TunitB TunitA
data2=np.copy(data)
for i in xrange(7):
    data2[:,i] = data[:,i+7]
    data2[:,i+7] = data[:, i]

data3 = np.concatenate((data, data2), axis=0)


# Pairs Within 1tunit
TSS_d=[]
CDS_group=[]
newData=[]
bed_d=[]
#Tunit1_size=[]
Q=[] # seperate the right Tunit into six regions. Q0: upstream of left TSS, Q1,Q2,Q3,Q4 of the left Tunit, Q5: downstream of PAS of left Tunit


for l in data3:
    l=list(l)
    #if l[4]!="S": # for OneS, Keep S on the right. Skip OneS with S on the left
    s1, e1 = l[1:3]
    s2, e2 = l[8:10]
    s1 = int(s1)
    s2 = int(s2)
    e1 = int(e1)
    e2 = int(e2)
    t1_size=abs(e1-s1)
    if (l[5]=="+" and min(abs(e2-s1),abs(e2-e1)) <= t1_size) or (l[5]=="-" and min(abs(s2-s1),abs(s2-e1)) <= t1_size):  # keep the pairs that TSS of right Tunit is within 1 Tunit distance of the left Tunit 
        Sbed1=(l[1], l[2], l[5])
        Sbed2=(l[8], l[9], l[12])
        #print l
        #print abs(s2-s1),abs(s2-e1), abs(e1-s1)
        #print "TSS distance", TSS_distance(Sbed1, Sbed2)
        newData.append(l)
        CDS_group.append(l[14]) #Con,Dis,OneS
        TSS_d.append(TSS_distance(Sbed1, Sbed2))
        #Tunit1_size.append(abs(e1-s1))
        if l[5]=="+":
            Q.append(round(float(e2 - s1)/t1_size, 3))
        else:
            Q.append(round(float(e1 - s2)/t1_size,3))


Q=np.array(Q)
CDS_group=np.array(CDS_group)
result={}

prefix='_'.join(infp.split("_")[0:2])
#print "CDS", "total","Q0","Q1", "Q2", "Q3","Q4","Q5"
q_cut=[-1.1, 0,0.25,0.5,0.75,1, 2.1]
q_name=["DivergentFromTSS", "Q1", "Q2", "Q3","Q4","PostPAS" ]
for g in ["Con","Dis", "OneS", "TwoS"]:
    q=Q[CDS_group==g]
    total=sum(CDS_group==g)*1.0
    for i in xrange(6):
    #if sum(CDS_group==g)>0 :
        print prefix, g, q_name[i], sum(np.logical_and(q>=q_cut[i], q<q_cut[i+1]))/max(total,1.0)
        
        #sum(np.logical_and(q>=0.25, q<0.5)),  sum(np.logical_and(q>=0.5, q<0.75)), \
        #sum(np.logical_and(q>=0.75, q<=1)),  sum(q>1)
    result[g]=[sum(CDS_group==g)*1.0, sum(q<0),  sum(np.logical_and(q>=0, q<0.25)), \
    sum(np.logical_and(q>=0.25, q<0.5)),  sum(np.logical_and(q>=0.5, q<0.75)), \
    sum(np.logical_and(q>=0.75, q<1)),  sum(q>=1)]


def barPlot(Con, Dis, pdfName, ylab="counts"):
    # data to plot
    n_groups = 6
    #Con =[0.35, 0.2, 0.02, 0.07, 0.05, 0.31]
    #Dis = [0.35, 0.2, 0.02, 0.07, 0.05, 0.31]
    
    # create plot
    fig, ax = plt.subplots()
    ax.set_title(pdfName)
    index = np.arange(n_groups)
    bar_width = 0.35
    opacity = 0.8
    
    rects1 = plt.bar(index, Con, bar_width,
    alpha=opacity,
    color='b',
    label='Con')
    
    rects2 = plt.bar(index + bar_width, Dis, bar_width,
    alpha=opacity,
    color='orange',
    label='Dis')
    
    plt.xlabel('TSS2 location at Tunit1')
    plt.ylabel(ylab)
    plt.title(pdfName)
    plt.xticks(index + bar_width/2, ('Before TSS', 'Q1', 'Q2', 'Q3', 'Q4', 'Post PAS'))
    plt.legend()
    
    plt.tight_layout()
    plt.savefig(pdfName)
    plt.show()
    


def barPlot3(Con, Dis, OneS, TwoS, pdfName, ylab="counts"):
    # data to plot
    n_groups = 6
    #Con =[0.35, 0.2, 0.02, 0.07, 0.05, 0.31]
    #Dis = [0.35, 0.2, 0.02, 0.07, 0.05, 0.31]
    
    # create plot
    fig, ax = plt.subplots()
    ax.set_title(pdfName)
    index = np.arange(n_groups)
    bar_width = 0.3
    opacity = 0.8
    
    rects1 = plt.bar(index, Con, bar_width,
    alpha=opacity,
    color='blue',
    label='Con')
    
    rects2 = plt.bar(index + bar_width, Dis, bar_width,
    alpha=opacity,
    color='orange',
    label='Dis')
    
    rects3 = plt.bar(index + 2*bar_width, OneS, bar_width,
    alpha=opacity,
    color='black',
    label='OneS')
    
    rects4 = plt.bar(index + 3*bar_width, TwoS, bar_width,
    alpha=opacity,
    color='red',
    label='TwoS')
    
    plt.xlabel('TSS2 location at Tunit1')
    plt.ylabel(ylab)
    #plt.title('Scores by person')
    plt.xticks(index + bar_width, ('Divergent from TSS', 'Q1', 'Q2', 'Q3', 'Q4', 'Post PAS'))
    plt.legend()
    
    plt.tight_layout()
    plt.savefig(pdfName)
    plt.show()

#pdfName="foo1.pdf"
barPlot(result["Con"][1:], result["Dis"][1:], pdfName+"_counts.pdf")
#barPlot(result["Con"][1:]/result["Con"][0], result["Dis"][1:]/result["Dis"][0], pdfName+"_fraction.pdf", ylab="fraction")

barPlot3(np.array(result["Con"][1:])/max(result["Con"][0],1), np.array(result["Dis"][1:])/max(1,result["Dis"][0]),np.array(result["OneS"][1:])/max(1,result["OneS"][0]), np.array(result["TwoS"][1:])/max(1,result["TwoS"][0]), pdfName+"_4fraction.pdf", ylab="fraction")


