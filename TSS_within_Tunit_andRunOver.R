#setwd("~/Box Sync/Danko_lab_work/F1_8Tissues/Tunit_pair_analysis/TSS_within_Tunit")
args=(commandArgs(TRUE))
setwd(args[1])
df_fp=args[2] #"TSS_within_Tunit_counts_report.txt"
pdf_name=args[3]#"TSS_within_Tunit_counts_report.pdf"


df=read.table(df_fp, header = F)
colnames(df)=c( "Tissue_cross", "CDS", "TSS_location","exp_fraction","control_fraction" )


df=df[grep("HT",df$Tissue_cross, invert = T),]

df=df[grep("ST",df$Tissue_cross, invert = T),]
df$TSS_location=factor(df$TSS_location, levels = c("-10K-5K" ,"-5K-TSS","TSS+5K", "5K-10K","10K-15K" ,"15K-20K", "20K-25K" ,"25K-30K" ))
df=df[!is.na(df$TSS_location),]
df$Tissue = substr(df$Tissue_cross,1,2)
df$fraction=df$exp_fraction
df$fraction=df$exp_fraction/df$control_fraction

#df$TSS_location[df$TSS_location=="Q0"] = "DivergentFromTSS"
#df$TSS_location[df$TSS_location=="Q5"] = "PostPAS"
head(df)

library(ggplot2)
library("RColorBrewer")
#display.brewer.all()
display.brewer.pal(n = 10, name = 'Paired')
paired_col=brewer.pal(n = 10, name = 'Paired')[c(1,2,7,8,3,4,5,6)]
single_col=paired_col[c(1,3,5)]

df$col="black"
df$col[df$CDS=="Con"]="#56B4E9"
df$col[df$CDS=="Dis"]="#E69F00"

col = brewer.pal(n = 8, name = 'Set1')
col = c(col,col)
col = sort(col) 
for (i in 1:length(levels(df$Tissue_cross))){
  df$col[df$Tissue_cross == levels(df$Tissue_cross)[i]]= col[i]
}


# grouped boxplot
pdf(pdf_name)
p=ggplot(df, aes(x=TSS_location, y=fraction, fill=CDS)) + 
  scale_color_brewer(palette="Paired") +
  #scale_fill_manual(values=single_col)+
  #scale_fill_manual(values=c("#56B4E9", "#E69F00","#999999", "red" )) +
  geom_boxplot() + 
  #geom_point(color=df$col)
  geom_jitter(shape=21, position=position_dodge(1))
p
dev.off()

pdf(pdf_name_2)
sub_df=df[df$CDS!="OneS",]
sub_df=df[grep("OneS",df$CDS, invert = T),]
subdf=sub_df[order(sub_df$TSS_location, sub_df$CDS),]
p=ggplot(sub_df, aes(x=TSS_location, y=fraction, fill=CDS)) + 
  #scale_color_brewer(palette="Dark2") +
  scale_fill_manual(values=c("white","yellow", "gray","brown"))+
  #scale_fill_manual(values=c("#56B4E9", "#E69F00","#999999", "red" )) +
  geom_boxplot() #+ 

p+scale_fill_manual(values=single_col)+ 
  geom_jitter(shape=21, position=position_dodge(1))
 # geom_jitter(shape=sub_df$Tissue_cross[order(sub_df$TSS_location, sub_df$CDS)], position=position_dodge(1))
p
p+ geom_jitter(aes(x=TSS_location, y=fraction, shape=Tissue, color =CDS)) + 
  scale_color_manual(values=c("blue","red","green","orange")) +
  scale_shape_manual(values=c(15,1,2,8,10,12,9,16,7))

    dev.off()


plot(df$TSS_location, df$fraction,)

pch=1:16














### old code below ###

# setwd("~/Box Sync/Danko_lab_work/F1_8Tissues/Tunit_pair_analysis/TSS_within_Tunit")
# df=read.table("TSS_within_Tunit_result.txt", header = T)
# 
# df$q0=df$Q0/df$total
# df$q1=df$Q1/df$total
# df$q2=df$Q2/df$total
# df$q3=df$Q3/df$total
# df$q4=df$Q4/df$total
# df$q5=df$Q5/df$total
# 
# # create a data frame
# gf=data.frame()
# Tissue_cross=rep(df$Tissue_cross, 6)
# CDS=rep(df$CDS,6)
# variety=rep(c("Q0","Q1","Q2","Q3","Q4","Q5"), each=length(df$Q0))
# fraction=c(df$q0,df$q1,df$q2,df$q3,df$q4,df$q5)
# gf=data.frame(Tissue_cross, CDS, variety,  fraction)
# 
# # examine
# df$q3[df$Tissue_cross=="ST_PB6" & df$CDS=="Dis"] == gf$fraction[gf$Tissue_cross=="ST_PB6" & gf$CDS=="Dis" & gf$variety=="Q3"]
# 
# 
# # library
# library(ggplot2)
# 
# # grouped boxplot
# p=ggplot(gf, aes(x=variety, y=fraction, fill=CDS)) + 
#   scale_fill_manual(values=c("#56B4E9", "#E69F00","#999999" )) +geom_boxplot()
# # p + geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(1), dotsize = 0.5)  
#  p+ geom_jitter(shape=21, position=position_dodge(1))
# 
# 
#  
# p
# p + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)
# #+  geom_dotplot(binaxis='y', stackdir='center',
# #                 position=position_dodge(1))

