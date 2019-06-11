#setwd("~/Box Sync/Danko_lab_work/F1_8Tissues/Tunit_pair_analysis/TSS_within_Tunit")
df=read.table("TSS_within_Tunit_counts_report.txt", header = F)
colnames(df)=c( "Tissue_cross", "CDS", "Q_name","fraction" )

df=df[grep("HT",df$Tissue_cross, invert = T),]
df$Q_name=factor(df$Q_name, levels = c("DivergentFromTSS","Q0", "Q1", "Q2", "Q3", "Q4", "Q5", "PostPAS"))
df$Q_name[df$Q_name=="Q0"] = "DivergentFromTSS"
df$Q_name[df$Q_name=="Q5"] = "PostPAS"


library(ggplot2)

# grouped boxplot
pdf("TSS_within_Tunit_counts_report.pdf")
p=ggplot(df, aes(x=Q_name, y=fraction, fill=CDS)) + 
  scale_fill_manual(values=c("#56B4E9", "#E69F00","#999999" )) +geom_boxplot()
p
p+ geom_jitter(shape=21, position=position_dodge(1))
dev.off()















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

