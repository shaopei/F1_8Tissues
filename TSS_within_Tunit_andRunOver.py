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


infp=argv[1] #"BN_MB6_paires_within_cluster1M.txt"
pdfName=argv[2]
#outfp=argv[2] #"BN_MB6_paires_within_cluster1M_overlapped_pairs.txt"
TSS_distance_max=50000
Tss_distance_bin=5000

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
data=[]
data2=[]

# Pairs Within 1tunit
TSS_d=[]
CDS_group=[]
newData=[]
bed_d=[]
#Tunit1_size=[]
Q=[] # seperate the right Tunit into six regions. Q0: upstream of left TSS, Q1,Q2,Q3,Q4 of the left Tunit, Q5: downstream of PAS of left Tunit
TSS_run_over=[]

newData_control=[]
TSS_run_over_control=[]
Q_control=[]
CDS_group_control=[]


for l in data3:
    l=list(l)
    s1, e1 = l[1:3]
    s2, e2 = l[8:10]
    s1 = int(s1)
    s2 = int(s2)
    e1 = int(e1)
    e2 = int(e2)
    t1_size=abs(e1-s1)
    t2_size=abs(e2-s2)
    if t1_size >= t2_size:
        Sbed1=(l[1], l[2], l[5])
        Sbed2=(l[8], l[9], l[12])
            #print l
            #print abs(s2-s1),abs(s2-e1), abs(e1-s1)
            #print "TSS distance", TSS_distance(Sbed1, Sbed2)
        tss_d=TSS_distance(Sbed1, Sbed2)
        if tss_d < TSS_distance_max:
            if l[4]!="S": #t1=M or P   and t2 = M/P/S (Con/Dis/OneS) 
                newData.append(l)
                CDS_group.append(l[14]) #Con,Dis,OneS
                TSS_d.append(tss_d)
                TSS_run_over.append(Sbed1_TSS_run_over_by_Sbed2(Sbed1, Sbed2))
                #Tunit1_size.append(abs(e1-s1))
                if l[5]=="+":
                    Q.append(round(float(e2 - s1)/Tss_distance_bin, 3))
                else:
                    Q.append(round(float(e1 - s2)/Tss_distance_bin,3))
            else: #t1=S, t2=M|P|S  
                    newData_control.append(l)
                    TSS_run_over_control.append(Sbed1_TSS_run_over_by_Sbed2(Sbed1, Sbed2))
                    CDS_group_control.append(l[14]) 
                    #Tunit1_size.append(abs(e1-s1))
                    if l[5]=="+":
                        Q_control.append(round(float(e2 - s1)/Tss_distance_bin, 3))
                    else:
                        Q_control.append(round(float(e1 - s2)/Tss_distance_bin,3))


Q=np.array(Q)
CDS_group=np.array(CDS_group)
TSS_run_over=np.array(TSS_run_over)
result={}

Q_control=np.array(Q_control)
TSS_run_over_control=np.array(TSS_run_over_control)
newData_control=np.array(newData_control)
CDS_group_control=np.array(CDS_group_control)

for g in ["Con","Dis", "OneS"]:
    result[g]=[]
    result[g+"RunOver"]=[]
    result[g+"NoRunOver"]=[]

prefix='_'.join(infp.split("_")[0:2])
#print "CDS", "total","Q0","Q1", "Q2", "Q3","Q4","Q5"
q_cut = range(-2,7) #[-1.1, 0,10.25,0.5,0.75,1, 2.1]
q_name= ["-10K-5K", "-5K-TSS", "TSS+5K", "5K-10K","10K-15K","15K-20K","20K-25K","25K-30K" ]

for i in xrange(len(q_name)):
    within_q_cut=np.logical_and( Q >= q_cut[i], Q < q_cut[i+1]) #within the range
    within_q_cut_control=np.logical_and( Q_control >= q_cut[i], Q_control < q_cut[i+1])
    total = sum(within_q_cut)*1.0  #t1=[M ||P]
    total_control = sum(within_q_cut_control)*1.0  # t1="S"
    #result[q_name[i]]=[]
    for g in ["Con","Dis", "OneS"]:
        #q = Q[np.logical_and(within_q_cut, CDS_group==g)]
        r=TSS_run_over[np.logical_and(within_q_cut, CDS_group==g)]
    #if sum(CDS_group==g)>0 :
        #print prefix, g+"_RunOver", q_name[i],  sum(r)/max(total,1.0)
        #print prefix, g+"_NoRunOver", q_name[i], (sum(np.logical_and(within_q_cut, CDS_group==g)) - sum(r))/max(total,1.0)
        
        #result[q_name[i]].append(sum(np.logical_and(within_q_cut, CDS_group==g)))
        result[g].append(round(sum(np.logical_and(within_q_cut, CDS_group==g))/max(total,1.0),5))
        result[g+"RunOver"].append(round((sum(r))/max(total,1.0),5))
        result[g+"NoRunOver"].append(round((sum(np.logical_and(within_q_cut, CDS_group==g)) - sum(r))/max(total,1.0),5))
        
    #                               P(t2="Con"|t1=[M ||P])                                       /     P(t2=[M||P]] | t1="S")
    #    print prefix, g, q_name[i], (sum(np.logical_and(within_q_cut, CDS_group==g)))/max(total,1.0), (sum(np.logical_and(within_q_cut_control, newData_control[:,11]!="S")))/max(total_control,1.0)
        
        #P(t1="Dis" |t2=[M||P])
        print prefix, g, q_name[i], sum(np.logical_and(within_q_cut, CDS_group==g)) ,sum(np.logical_and(within_q_cut, CDS_group!="OneS")) + sum(np.logical_and(within_q_cut_control, CDS_group_control =="OneS"))
        #P(t1="Dis" |t2=[M||P] && t2 cross t1 TSS)
        print prefix, g+"RunOver", q_name[i], sum(np.logical_and(np.logical_and(within_q_cut, CDS_group==g), TSS_run_over)) ,sum(np.logical_and(np.logical_and(within_q_cut, CDS_group!="OneS"), TSS_run_over))+ sum(np.logical_and(np.logical_and(within_q_cut_control, CDS_group_control =="OneS"), TSS_run_over_control))



def barPlot(Con, Dis, pdfName, ylab="frequency"):
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
    plt.xticks(index + bar_width/2, ("-10K-5K", "-5K-TSS", "TSS+5K", "5K-10K","10K-15K","15K-20K"))
    plt.legend()
    
    plt.tight_layout()
    plt.savefig(pdfName)
    plt.show()
    


def barPlotRunOver(ConNoRunOver, ConRunOver, DisNoRunOver, DisRunOver, pdfName, ylab="fraction"):
    # data to plot
    n_groups = 6
    #Con =[0.35, 0.2, 0.02, 0.07, 0.05, 0.31]
    #Dis = [0.35, 0.2, 0.02, 0.07, 0.05, 0.31]
    
    # create plot
    fig, ax = plt.subplots()
    ax.set_title(pdfName)
    index = np.arange(n_groups)
    bar_width = 0.15
    opacity = 0.8
    
    rects1 = plt.bar(index, ConNoRunOver, bar_width,
    alpha=opacity,
    color='b',
    label='Con NoRunOver')
    
    rects2 = plt.bar(index + bar_width, ConRunOver, bar_width,
    alpha=opacity,
    color='green',
    label='Con RunOver')
    
    rects3 = plt.bar(index + 2*bar_width, DisNoRunOver, bar_width,
    alpha=opacity,
    color='orange',
    label='Dis NoRunOver')
    
    rects4 = plt.bar(index + 3*bar_width, DisRunOver, bar_width,
    alpha=opacity,
    color='red',
    label='Dis RunOver')
    
    plt.xlabel('TSS2 location at Tunit1')
    plt.ylabel(ylab)
    plt.title(pdfName)
    plt.xticks(index + bar_width, ("-10K-5K", "-5K-TSS", "TSS+5K", "5K-10K","10K-15K","15K-20K"))
    plt.legend()
    
    plt.tight_layout()
    plt.savefig(pdfName)
    plt.show()
    


def barPlot3(Con, Dis, OneS, pdfName, ylab="counts"):
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
    
    plt.xlabel('TSS2 location at Tunit1')
    plt.ylabel(ylab)
    #plt.title('Scores by person')
    plt.xticks(index + bar_width,  ("-10K-5K", "-5K-TSS", "TSS+5K", "5K-10K","10K-15K","15K-20K"))
    plt.legend()
    
    plt.tight_layout()
    plt.savefig(pdfName)
    plt.show()

#pdfName="foo1.pdf"
#barPlot(result["Con"], result["Dis"], pdfName+"_frequency.pdf")
#barPlot3(result["Con"], result["Dis"], result["OneS"], pdfName+"_3frequency.pdf")
#barPlot3(np.array(result["Con"][1:])/max(result["Con"][0],1), np.array(result["Dis"][1:])/max(1,result["Dis"][0]),np.array(result["OneS"][1:])/max(1,result["OneS"][0]), pdfName+"_fraction.pdf", ylab="fraction")
#barPlotRunOver(result["ConNoRunOver"], result["ConRunOver"], result["DisNoRunOver"], result["DisRunOver"], pdfName, ylab="fraction")

