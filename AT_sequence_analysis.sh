cd /workdir/sc2457/F1_Tissues/PolyA_Allele-specific-manuscript

# get the M,P identify fi AlleleHMMblocks to B6 and CAST

Organ=LV
# from MB6
intersectBed -wo -s -f 1 -r -a ${Organ}_AT_4tunitIntersectNativeHMM_AlleleHMM.bed -b HMM_bed/${Organ}_MB6_HMM_*.bed \
| awk 'BEGIN {OFS="\t"}($4==$17){print $1,$2,$3,$4,$5,$6, $7, $8, $9, $10, $11, $12}' \
| awk 'BEGIN {OFS="\t"} (substr($4,1,1)=="M")  {print $0, "B6"} (substr($4,1,1)=="P")  {print $0, "CAST"}' > ${Organ}_AT_4tunitIntersectNativeHMM_AlleleHMM_strain.bed

# from PB6
intersectBed -wo -s -f 1 -r -a ${Organ}_AT_4tunitIntersectNativeHMM_AlleleHMM.bed -b HMM_bed/${Organ}_PB6_HMM_*.bed \
| awk 'BEGIN {OFS="\t"}($4==$17){print $1,$2,$3,$4,$5,$6, $7, $8, $9, $10, $11, $12}' \
| awk 'BEGIN {OFS="\t"} (substr($4,1,1)=="P")  {print $0, "B6"} (substr($4,1,1)=="M")  {print $0, "CAST"}' >> ${Organ}_AT_4tunitIntersectNativeHMM_AlleleHMM_strain.bed

# get the intersectRegion
intersectBed -wb -s -a <(cat ${Organ}_AT_4tunitIntersectNativeHMM_AlleleHMM_strain.bed | cut -f 7-12) -b ${Organ}_AT_4tunitIntersectNativeHMM_AlleleHMM_strain.bed |sort |uniq >${Organ}_AT_4tunitIntersectNativeHMM_intersectRegion_strain.bed

# get the 1st bp of AlleleHMM blocks inside Tunits
cat ${Organ}_AT_4tunitIntersectNativeHMM_intersectRegion_strain.bed | awk 'BEGIN {OFS="\t"} ($6=="+") {print $1, $2, $2+1, $4, $19, $6} ($6=="-") {print $1, $3-1, $3, $4, $19, $6}' > ${Organ}_AT_4tunitIntersectNativeHMM_intersectRegion_strain_temp1.bed


ln -s /workdir/sc2457/F1_Tissues/TSN_SingleBaseRunOn/P.CAST.EiJ_M.C57BL.6J_*aternal_all.fa* .

j=${Organ}_AT_4tunitIntersectNativeHMM_intersectRegion_strain_temp1
d=100
#if [ ! -f ${j}_+-${d}_High_LowAlleleSeq.bed ]; then
 # get sequence from maternal genome
 bedtools getfasta -s -fi P.CAST.EiJ_M.C57BL.6J_maternal_all.fa -bed <(cat ${j}.bed |awk -v d=$d  '{OFS="\t";p="_maternal"} {print substr($1,4)p, $2-d, $3+d, $4,$5,$6}')  | grep -v \> > ${j}_P.CAST.EiJ_M.C57BL.6J_maternal.txt &
 # get sequence from paternal genome
 bedtools getfasta -s -fi P.CAST.EiJ_M.C57BL.6J_paternal_all.fa -bed <(cat ${j}.bed |awk -v d=$d  '{OFS="\t";p="_paternal"} {print substr($1,4)p, $2-d, $3+d, $4,$5,$6}')  | grep -v \> > ${j}_P.CAST.EiJ_M.C57BL.6J_paternal.txt &
 wait
 paste ${j}.bed  ${j}_P.CAST.EiJ_M.C57BL.6J_maternal.txt ${j}_P.CAST.EiJ_M.C57BL.6J_paternal.txt > ${j}_+-${d}_B6_CAST_Seq.bed 
 cat ${j}_+-${d}_B6_CAST_Seq.bed  | awk '{OFS="\t"} ($5=="B6") {print $1,$2,$3,$4,$5, $6,"." ,$7, $8} 
 ($5=="CAST") {print  $1,$2,$3,$4,$5, $6,"." ,$8, $7} ' > ${j}_+-${d}_Long_ShortAlleleSeq.bed 

ln -s /workdir/sc2457/mouse_AlleleSpecific/mouse_genome.sanger.ac.uk/REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/P.CAST_M.B6_indelsNsnps_CAST.bam.alleleDBInput.snp.sorted.bed .
Organ=BN
# identify 1SNP identity
intersectBed -wo -a ${Organ}_AT_4tunitIntersectNativeHMM_intersectRegion_strain_temp1.bed -b P.CAST_M.B6_indelsNsnps_CAST.bam.alleleDBInput.snp.sorted.bed > ${Organ}_AT_4tunitIntersectNativeHMM_intersectRegion_strain_temp1_SNP.bed

# identify 2nd SNP downstream of 1st bp of AT window
bedtools closest -io -iu -D a -a <(cat ${Organ}_AT_4tunitIntersectNativeHMM_intersectRegion_strain_temp1_SNP.bed| sort -k1,1 -k2,2n) -b <(cat P.CAST_M.B6_indelsNsnps_CAST.bam.alleleDBInput.snp.sorted.bed| sort -k1,1 -k2,2n) | uniq > ${Organ}_AT_4tunitIntersectNativeHMM_intersectRegion_strain_temp1_2SNP.bed

# identify ATs with 1st bp short allele C, long allele A,T,G
# with 2nd SNPs within 20bp
cat ${Organ}_AT_4tunitIntersectNativeHMM_intersectRegion_strain_temp1_2SNP.bed \
| awk 'BEGIN {OFS="\t"}($5=="B6" && $6=="+" && substr($10,2,2)=="C") {print $0} ($5=="B6" && $6=="-" && substr($10,2,2)=="G") {print $0} 
($5=="CAST" && $6=="+" && substr($10,1,1)=="C") {print $0} ($5=="CAST" && $6=="-" && substr($10,1,1)=="G") {print $0}' \
| awk 'BEGIN {OFS="\t";c="_SNP2"}($16+1 > 0 && $16+0 <20) {print $1, $13, $14, $4c, $5, $6 ,$10, $15}'  \
>  ${Organ}_AT_4tunitIntersectNativeHMM_intersectRegion_strain_temp1_SNP2.bed

j=${Organ}_AT_4tunitIntersectNativeHMM_intersectRegion_strain_temp1_SNP2
d=200
#if [ ! -f ${j}_+-${d}_High_LowAlleleSeq.bed ]; then
 # get sequence from maternal genome
 bedtools getfasta -s -fi P.CAST.EiJ_M.C57BL.6J_maternal_all.fa -bed <(cat ${j}.bed |awk -v d=$d  '{OFS="\t";p="_maternal"} {print substr($1,4)p, $2-d, $3+d, $4,$5,$6}')  | grep -v \> > ${j}_P.CAST.EiJ_M.C57BL.6J_maternal.txt &
 # get sequence from paternal genome
 bedtools getfasta -s -fi P.CAST.EiJ_M.C57BL.6J_paternal_all.fa -bed <(cat ${j}.bed |awk -v d=$d  '{OFS="\t";p="_paternal"} {print substr($1,4)p, $2-d, $3+d, $4,$5,$6}')  | grep -v \> > ${j}_P.CAST.EiJ_M.C57BL.6J_paternal.txt &
 wait
 paste ${j}.bed  ${j}_P.CAST.EiJ_M.C57BL.6J_maternal.txt ${j}_P.CAST.EiJ_M.C57BL.6J_paternal.txt > ${j}_+-${d}_B6_CAST_Seq.bed 
 cat ${j}_+-${d}_B6_CAST_Seq.bed  | awk '{OFS="\t"} ($5=="B6") {print $1,$2,$3,$4,$5, $6,$7, $8, $9, $10} 
 ($5=="CAST") {print  $1,$2,$3,$4,$5, $6, $7 ,$8, $10, $9} ' > ${j}_+-${d}_Long_ShortAlleleSeq.bed 


# check all 2nd SNPs  
Organ=BN
cat ${Organ}_AT_4tunitIntersectNativeHMM_intersectRegion_strain_temp1_2SNP.bed \
| awk 'BEGIN {OFS="\t";c="_SNP2"}($16+1 > 0) {print $1, $13, $14, $4c, $5, $6 ,$10, $15}'  \
>  ${Organ}_AT_4tunitIntersectNativeHMM_intersectRegion_strain_temp1_SNP2.bed

j=${Organ}_AT_4tunitIntersectNativeHMM_intersectRegion_strain_temp1_SNP2
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







