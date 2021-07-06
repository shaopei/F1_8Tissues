cd /workdir/sc2457/F1_Tissues/PolyA_Allele-specific-manuscript

# get the M,P identify fi AlleleHMMblocks to B6 and CAST

Organ=LV
# annotate which stran the AlleleHMM blocks biased to
# from MB6
# -f 1 -r require complete overlap of -a and -b
intersectBed -wo -s -f 1 -r -a ${Organ}_AT_4tunitIntersectNativeHMM_AlleleHMM.bed -b HMM_bed/${Organ}_MB6_HMM_*.bed \
| awk 'BEGIN {OFS="\t"}($4==$17){print $1,$2,$3,$4,$5,$6, $7, $8, $9, $10, $11, $12}' \
| awk 'BEGIN {OFS="\t"} (substr($4,1,1)=="M")  {print $0, "B6"} (substr($4,1,1)=="P")  {print $0, "CAST"}' | sort-bed - > ${Organ}_AT_4tunitIntersectNativeHMM_AlleleHMM_strain.temp.bed

# from PB6
intersectBed -wo -s -f 1 -r -a ${Organ}_AT_4tunitIntersectNativeHMM_AlleleHMM.bed -b HMM_bed/${Organ}_PB6_HMM_*.bed \
| awk 'BEGIN {OFS="\t"}($4==$17){print $1,$2,$3,$4,$5,$6, $7, $8, $9, $10, $11, $12}' \
| awk 'BEGIN {OFS="\t"} (substr($4,1,1)=="P")  {print $0, "B6"} (substr($4,1,1)=="M")  {print $0, "CAST"}' | sort-bed - >> ${Organ}_AT_4tunitIntersectNativeHMM_AlleleHMM_strain.temp.bed

cat ${Organ}_AT_4tunitIntersectNativeHMM_AlleleHMM_strain.temp.bed |sort-bed - |uniq > ${Organ}_AT_4tunitIntersectNativeHMM_AlleleHMM_strain.bed  

# get the intersectRegion
# use ($4==$16)  to make sure the intersect coms from correct pairs
intersectBed -wb -s -a <(cat ${Organ}_AT_4tunitIntersectNativeHMM_AlleleHMM_strain.bed | cut -f 7-12) -b ${Organ}_AT_4tunitIntersectNativeHMM_AlleleHMM_strain.bed | awk 'BEGIN {OFS="\t"} ($4==$16) {print $0}' | sort-bed - |uniq >${Organ}_AT_4tunitIntersectNativeHMM_intersectRegion_strain.bed

# get the 1st bp of AlleleHMM blocks inside Tunits
# 7-12 is the location of AT
cat ${Organ}_AT_4tunitIntersectNativeHMM_intersectRegion_strain.bed | awk 'BEGIN {OFS="\t"} ($6=="+") {print $1, $2, $2+1, $4, $19, $6, $1, $2, $3, $4, $5, $6} ($6=="-") {print $1, $3-1, $3, $4, $19, $6, $1, $2, $3, $4, $5, $6}' |sort-bed - > ${Organ}_AT_4tunitIntersectNativeHMM_intersectRegion_strain_1stBp.bed


ln -s /workdir/sc2457/F1_Tissues/TSN_SingleBaseRunOn/P.CAST.EiJ_M.C57BL.6J_*aternal_all.fa* .

j=${Organ}_AT_4tunitIntersectNativeHMM_intersectRegion_strain_1stBp
d=30
#if [ ! -f ${j}_+-${d}_High_LowAlleleSeq.bed ]; then
 # get sequence from maternal genome
 bedtools getfasta -s -fi P.CAST.EiJ_M.C57BL.6J_maternal_all.fa -bed <(cat ${j}.bed |awk -v d=$d  '{OFS="\t";p="_maternal"} {print substr($1,4)p, $2-d, $3+d, $4,$5,$6}')  | grep -v \> > ${j}_P.CAST.EiJ_M.C57BL.6J_maternal.txt &
 # get sequence from paternal genome
 bedtools getfasta -s -fi P.CAST.EiJ_M.C57BL.6J_paternal_all.fa -bed <(cat ${j}.bed |awk -v d=$d  '{OFS="\t";p="_paternal"} {print substr($1,4)p, $2-d, $3+d, $4,$5,$6}')  | grep -v \> > ${j}_P.CAST.EiJ_M.C57BL.6J_paternal.txt &
 wait
 paste ${j}.bed  ${j}_P.CAST.EiJ_M.C57BL.6J_maternal.txt ${j}_P.CAST.EiJ_M.C57BL.6J_paternal.txt > ${j}_+-${d}_B6_CAST_Seq.bed 
 cat ${j}_+-${d}_B6_CAST_Seq.bed  | awk '{OFS="\t"} ($5=="B6") {print $1,$2,$3,$4,$5, $6,"." ,$13, $14} 
 ($5=="CAST") {print  $1,$2,$3,$4,$5, $6,"." ,$14, $13} ' > ${j}_+-${d}_Long_ShortAlleleSeq.bed 
# $5 indicates which strain is long allele

ln -s /workdir/sc2457/mouse_AlleleSpecific/mouse_genome.sanger.ac.uk/REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/P.CAST_M.B6_indelsNsnps_CAST.bam.alleleDBInput.snp.sorted.bed .
cat P.CAST_M.B6_indelsNsnps_CAST.bam.alleleDBInput.snp.sorted.bed | sort-bed - > P.CAST_M.B6_indelsNsnps_CAST.bam.alleleDBInput.snp.sort-bed.bed
Organ=LV
## identify 1SNP identity
intersectBed -wo -a ${Organ}_AT_4tunitIntersectNativeHMM_intersectRegion_strain_1stBp.bed -b P.CAST_M.B6_indelsNsnps_CAST.bam.alleleDBInput.snp.sort-bed.bed > ${Organ}_AT_4tunitIntersectNativeHMM_intersectRegion_strain_1stBp_SNP.bed
# 1-6 1st bp
# 7-12 AT
# 13-16 1st SNP identity
# 17 region overlap =1

## identify 2nd SNP downstream of 1st bp of AT window
bedtools closest -io -iu -D a -a <(cat ${Organ}_AT_4tunitIntersectNativeHMM_intersectRegion_strain_1stBp_SNP.bed| sort-bed -) -b P.CAST_M.B6_indelsNsnps_CAST.bam.alleleDBInput.snp.sort-bed.bed  > ${Organ}_AT_4tunitIntersectNativeHMM_intersectRegion_strain_1stBp_2ndSNP.bed
# 1-6 1st bp
# 7-12 AT
# 13-16 SNP1 identity
# 17 region overlap =1
# 18-21 SNP2 identity
# 22 distance(SNP1, SNP2)

## identify 3rd SNP downstream of 1st bp of AT window
bedtools closest -io -iu -D a -a <(cat ${Organ}_AT_4tunitIntersectNativeHMM_intersectRegion_strain_1stBp_2ndSNP.bed| awk 'BEGIN {OFS="\t"} {print $18, $19, $20, $4, $5, $6, $7, $8, $9, $10, $5, $12, $13, $14 ,$15, $16, $17, $18, $19, $20, $21, $22}'| sort-bed -) -b P.CAST_M.B6_indelsNsnps_CAST.bam.alleleDBInput.snp.sort-bed.bed \
 > ${Organ}_AT_4tunitIntersectNativeHMM_intersectRegion_strain_1stBp_2ndSNP_3rdSNP.bed
# 1-6 2nd SNP
# 7-12 AT
# 13-16 SNP1 identity
# 17 region overlap =1
# 18-21 SNP2 identity
# 22 distance(SNP1, SNP2)
# 23-26 SNP2 identity
# 27 distance(SNP2, SNP3)
## identify 4th SNP downstream of 1st bp of AT window
bedtools closest -io -iu -D a -a <(cat ${Organ}_AT_4tunitIntersectNativeHMM_intersectRegion_strain_1stBp_2ndSNP_3rdSNP.bed| awk 'BEGIN {OFS="\t"} {print $23, $24, $25, $4, $5, $6, $0}' | cut -f 1-6,13-| sort-bed -) -b P.CAST_M.B6_indelsNsnps_CAST.bam.alleleDBInput.snp.sort-bed.bed \
  > ${Organ}_AT_4tunitIntersectNativeHMM_intersectRegion_strain_1stBp_2ndSNP_3rdSNP_4thSNP.bed
# 1-6 3rd SNP
# 28-31 SNP3 identity
# 32 distance(SNP3, SNP4)

## identify 5th SNP downstream of 1st bp of AT window
bedtools closest -io -iu -D a -a <(cat ${Organ}_AT_4tunitIntersectNativeHMM_intersectRegion_strain_1stBp_2ndSNP_3rdSNP_4thSNP.bed| awk 'BEGIN {OFS="\t"} {print $28, $29, $30, $4, $5, $6, $0}' | cut -f 1-6,13-| sort-bed -) -b P.CAST_M.B6_indelsNsnps_CAST.bam.alleleDBInput.snp.sort-bed.bed  \
> ${Organ}_AT_4tunitIntersectNativeHMM_intersectRegion_strain_1stBp_2ndSNP_3rdSNP_4thSNP_5thSNP.bed
#

# check all 2nd SNPs  
Organ=LV
cat ${Organ}_AT_4tunitIntersectNativeHMM_intersectRegion_strain_1stBp_2ndSNP_3rdSNP.bed \
| awk 'BEGIN {OFS="\t";c="_SNP2"} {print $1, $2, $3, $4c, $5, $6 }'  \
> ${Organ}_AT_4tunitIntersectNativeHMM_intersectRegion_strain_1stBp_SNP2.bed 
#>  ${Organ}_AT_4tunitIntersectNativeHMM_intersectRegion_strain_1stBp_SNP2.bed

j=${Organ}_AT_4tunitIntersectNativeHMM_intersectRegion_strain_1stBp_SNP2
d=100
#if [ ! -f ${j}_+-${d}_High_LowAlleleSeq.bed ]; then
 # get sequence from maternal genome
 bedtools getfasta -s -fi P.CAST.EiJ_M.C57BL.6J_maternal_all.fa -bed <(cat ${j}.bed|awk -v d=$d  '{OFS="\t";p="_maternal"} {print substr($1,4)p, $2-d, $3+d, $4,$5,$6}')  | grep -v \> > ${j}_P.CAST.EiJ_M.C57BL.6J_maternal.txt &
 # get sequence from paternal genome
 bedtools getfasta -s -fi P.CAST.EiJ_M.C57BL.6J_paternal_all.fa -bed <(cat ${j}.bed|awk -v d=$d  '{OFS="\t";p="_paternal"} {print substr($1,4)p, $2-d, $3+d, $4,$5,$6}')  | grep -v \> > ${j}_P.CAST.EiJ_M.C57BL.6J_paternal.txt &
 wait
 paste ${j}.bed  ${j}_P.CAST.EiJ_M.C57BL.6J_maternal.txt ${j}_P.CAST.EiJ_M.C57BL.6J_paternal.txt > ${j}_+-${d}_B6_CAST_Seq.bed 
 cat ${j}_+-${d}_B6_CAST_Seq.bed  | awk '{OFS="\t"} ($5=="B6") {print $1,$2,$3,$4,$5, $6,$7, $8} 
 ($5=="CAST") {print  $1,$2,$3,$4,$5, $6, $8 ,$7} ' > ${j}_+-${d}_Long_ShortAlleleSeq.bed 


#HERE

# identify ATs with 1st bp short allele C, long allele A,T,G
cat ${Organ}_AT_4tunitIntersectNativeHMM_intersectRegion_strain_1stBp_2ndSNP_3rdSNP_4thSNP_5thSNP.bed \
| cut -f 7- \
| awk 'BEGIN {OFS="\t"}($5=="B6" && $6=="+" && substr($10,2,2)=="C") {print $0} ($5=="B6" && $6=="-" && substr($10,2,2)=="G") {print $0} 
($5=="CAST" && $6=="+" && substr($10,1,1)=="C") {print $0} ($5=="CAST" && $6=="-" && substr($10,1,1)=="G") {print $0}' \
> ${Organ}_AT_4tunitIntersectNativeHMM_intersectRegion_strain_1stBpShortC_2345thSNP.bed

# 1st SNP
j=${Organ}_AT_4tunitIntersectNativeHMM_intersectRegion_strain_1stBpShortC_SNP1
cat ${Organ}_AT_4tunitIntersectNativeHMM_intersectRegion_strain_1stBpShortC_2345thSNP.bed \
| awk 'BEGIN {OFS="\t";c="_SNP1"} {print $7, $8, $9, $4c, $5, $6 ,$10, $10}'  \
>  ${j}.bed 
d=100
#if [ ! -f ${j}_+-${d}_High_LowAlleleSeq.bed ]; then
 # get sequence from maternal genome
 bedtools getfasta -s -fi P.CAST.EiJ_M.C57BL.6J_maternal_all.fa -bed <(cat ${j}.bed |awk -v d=$d  '{OFS="\t";p="_maternal"} {print substr($1,4)p, $2-d, $3+d, $4,$5,$6}')  | grep -v \> > ${j}_P.CAST.EiJ_M.C57BL.6J_maternal.txt &
 # get sequence from paternal genome
 bedtools getfasta -s -fi P.CAST.EiJ_M.C57BL.6J_paternal_all.fa -bed <(cat ${j}.bed |awk -v d=$d  '{OFS="\t";p="_paternal"} {print substr($1,4)p, $2-d, $3+d, $4,$5,$6}')  | grep -v \> > ${j}_P.CAST.EiJ_M.C57BL.6J_paternal.txt &
 wait
 paste ${j}.bed  ${j}_P.CAST.EiJ_M.C57BL.6J_maternal.txt ${j}_P.CAST.EiJ_M.C57BL.6J_paternal.txt > ${j}_+-${d}_B6_CAST_Seq.bed 
 cat ${j}_+-${d}_B6_CAST_Seq.bed  | awk '{OFS="\t"} ($5=="B6") {print $1,$2,$3,$4,$5, $6,$7, $8, $9, $10} 
 ($5=="CAST") {print  $1,$2,$3,$4,$5, $6, $7 ,$8, $10, $9} ' > ${j}_+-${d}_Long_ShortAlleleSeq.bed 



#2nd SNP
cat ${Organ}_AT_4tunitIntersectNativeHMM_intersectRegion_strain_1stBpShortC_2345thSNP.bed \
| awk 'BEGIN {OFS="\t";c="_SNP2"}{print $12, $13, $14, $4c, $5, $6 ,$10, $15}'  \
>   ${Organ}_AT_4tunitIntersectNativeHMM_intersectRegion_strain_1stBpShortC_SNP2.bed

j=${Organ}_AT_4tunitIntersectNativeHMM_intersectRegion_strain_1stBpShortC_SNP2
d=100
#if [ ! -f ${j}_+-${d}_High_LowAlleleSeq.bed ]; then
 # get sequence from maternal genome
 bedtools getfasta -s -fi P.CAST.EiJ_M.C57BL.6J_maternal_all.fa -bed <(cat ${j}.bed |awk -v d=$d  '{OFS="\t";p="_maternal"} {print substr($1,4)p, $2-d, $3+d, $4,$5,$6}')  | grep -v \> > ${j}_P.CAST.EiJ_M.C57BL.6J_maternal.txt &
 # get sequence from paternal genome
 bedtools getfasta -s -fi P.CAST.EiJ_M.C57BL.6J_paternal_all.fa -bed <(cat ${j}.bed |awk -v d=$d  '{OFS="\t";p="_paternal"} {print substr($1,4)p, $2-d, $3+d, $4,$5,$6}')  | grep -v \> > ${j}_P.CAST.EiJ_M.C57BL.6J_paternal.txt &
 wait
 paste ${j}.bed  ${j}_P.CAST.EiJ_M.C57BL.6J_maternal.txt ${j}_P.CAST.EiJ_M.C57BL.6J_paternal.txt > ${j}_+-${d}_B6_CAST_Seq.bed 
 cat ${j}_+-${d}_B6_CAST_Seq.bed  | awk '{OFS="\t"} ($5=="B6") {print $1,$2,$3,$4,$5, $6,$7, $8, $9, $10} 
 ($5=="CAST") {print  $1,$2,$3,$4,$5, $6, $7 ,$8, $10, $9} ' > ${j}_+-${d}_Long_ShortAlleleSeq.bed 

#3rd SNP
j=${Organ}_AT_4tunitIntersectNativeHMM_intersectRegion_strain_1stBpShortC_SNP3
cat ${Organ}_AT_4tunitIntersectNativeHMM_intersectRegion_strain_1stBpShortC_2345thSNP.bed \
| awk 'BEGIN {OFS="\t";c="_SNP3"} {print $17, $18, $19, $4c, $5, $6 ,$10, $20}'  \
>   ${j}.bed

d=100
#if [ ! -f ${j}_+-${d}_High_LowAlleleSeq.bed ]; then
 # get sequence from maternal genome
 bedtools getfasta -s -fi P.CAST.EiJ_M.C57BL.6J_maternal_all.fa -bed <(cat ${j}.bed |awk -v d=$d  '{OFS="\t";p="_maternal"} {print substr($1,4)p, $2-d, $3+d, $4,$5,$6}')  | grep -v \> > ${j}_P.CAST.EiJ_M.C57BL.6J_maternal.txt &
 # get sequence from paternal genome
 bedtools getfasta -s -fi P.CAST.EiJ_M.C57BL.6J_paternal_all.fa -bed <(cat ${j}.bed |awk -v d=$d  '{OFS="\t";p="_paternal"} {print substr($1,4)p, $2-d, $3+d, $4,$5,$6}')  | grep -v \> > ${j}_P.CAST.EiJ_M.C57BL.6J_paternal.txt &
 wait
 paste ${j}.bed  ${j}_P.CAST.EiJ_M.C57BL.6J_maternal.txt ${j}_P.CAST.EiJ_M.C57BL.6J_paternal.txt > ${j}_+-${d}_B6_CAST_Seq.bed 
 cat ${j}_+-${d}_B6_CAST_Seq.bed  | awk '{OFS="\t"} ($5=="B6") {print $1,$2,$3,$4,$5, $6,$7, $8, $9, $10} 
 ($5=="CAST") {print  $1,$2,$3,$4,$5, $6, $7 ,$8, $10, $9} ' > ${j}_+-${d}_Long_ShortAlleleSeq.bed 

#4th SNP
j=${Organ}_AT_4tunitIntersectNativeHMM_intersectRegion_strain_1stBpShortC_SNP4
cat ${Organ}_AT_4tunitIntersectNativeHMM_intersectRegion_strain_1stBpShortC_2345thSNP.bed \
| awk 'BEGIN {OFS="\t";c="_SNP4"} {print $22, $23, $24, $4c, $5, $6 ,$10, $25}'  \
>   ${j}.bed

d=100
#if [ ! -f ${j}_+-${d}_High_LowAlleleSeq.bed ]; then
 # get sequence from maternal genome
 bedtools getfasta -s -fi P.CAST.EiJ_M.C57BL.6J_maternal_all.fa -bed <(cat ${j}.bed |awk -v d=$d  '{OFS="\t";p="_maternal"} {print substr($1,4)p, $2-d, $3+d, $4,$5,$6}')  | grep -v \> > ${j}_P.CAST.EiJ_M.C57BL.6J_maternal.txt &
 # get sequence from paternal genome
 bedtools getfasta -s -fi P.CAST.EiJ_M.C57BL.6J_paternal_all.fa -bed <(cat ${j}.bed |awk -v d=$d  '{OFS="\t";p="_paternal"} {print substr($1,4)p, $2-d, $3+d, $4,$5,$6}')  | grep -v \> > ${j}_P.CAST.EiJ_M.C57BL.6J_paternal.txt &
 wait
 paste ${j}.bed  ${j}_P.CAST.EiJ_M.C57BL.6J_maternal.txt ${j}_P.CAST.EiJ_M.C57BL.6J_paternal.txt > ${j}_+-${d}_B6_CAST_Seq.bed 
 cat ${j}_+-${d}_B6_CAST_Seq.bed  | awk '{OFS="\t"} ($5=="B6") {print $1,$2,$3,$4,$5, $6,$7, $8, $9, $10} 
 ($5=="CAST") {print  $1,$2,$3,$4,$5, $6, $7 ,$8, $10, $9} ' > ${j}_+-${d}_Long_ShortAlleleSeq.bed 

#5th SNP
j=${Organ}_AT_4tunitIntersectNativeHMM_intersectRegion_strain_1stBpShortC_SNP5
cat ${Organ}_AT_4tunitIntersectNativeHMM_intersectRegion_strain_1stBpShortC_2345thSNP.bed \
| awk 'BEGIN {OFS="\t";c="_SNP5"} {print $27, $28, $29, $4c, $5, $6 ,$10, $30}'  \
>   ${j}.bed

d=100
#if [ ! -f ${j}_+-${d}_High_LowAlleleSeq.bed ]; then
 # get sequence from maternal genome
 bedtools getfasta -s -fi P.CAST.EiJ_M.C57BL.6J_maternal_all.fa -bed <(cat ${j}.bed |awk -v d=$d  '{OFS="\t";p="_maternal"} {print substr($1,4)p, $2-d, $3+d, $4,$5,$6}')  | grep -v \> > ${j}_P.CAST.EiJ_M.C57BL.6J_maternal.txt &
 # get sequence from paternal genome
 bedtools getfasta -s -fi P.CAST.EiJ_M.C57BL.6J_paternal_all.fa -bed <(cat ${j}.bed |awk -v d=$d  '{OFS="\t";p="_paternal"} {print substr($1,4)p, $2-d, $3+d, $4,$5,$6}')  | grep -v \> > ${j}_P.CAST.EiJ_M.C57BL.6J_paternal.txt &
 wait
 paste ${j}.bed  ${j}_P.CAST.EiJ_M.C57BL.6J_maternal.txt ${j}_P.CAST.EiJ_M.C57BL.6J_paternal.txt > ${j}_+-${d}_B6_CAST_Seq.bed 
 cat ${j}_+-${d}_B6_CAST_Seq.bed  | awk '{OFS="\t"} ($5=="B6") {print $1,$2,$3,$4,$5, $6,$7, $8, $9, $10} 
 ($5=="CAST") {print  $1,$2,$3,$4,$5, $6, $7 ,$8, $10, $9} ' > ${j}_+-${d}_Long_ShortAlleleSeq.bed 





