cd /workdir/sc2457/F1_Tissues/PolyA_Allele-specific-manuscript

# get the M,P identify fi AlleleHMMblocks to B6 and CAST

Organ=BN
# from MB6
intersectBed -wo -s -f 1 -r -a ${Organ}_AT_4tunitIntersectNativeHMM_AlleleHMM.bed -b HMM_bed/${Organ}_MB6_HMM_*.bed \
| awk 'BEGIN {OFS="\t"}($4==$17){print $1,$2,$3,$4,$5,$6, $7, $8, $9, $10, $11, $12}' \
| awk 'BEGIN {OFS="\t"} (substr($4,1,1)=="M")  {print $0, "B6"} (substr($4,1,1)=="P")  {print $0, "CAST"}' > ${Organ}_AT_4tunitIntersectNativeHMM_AlleleHMM_strain.bed

# from PB6
intersectBed -wo -s -f 1 -r -a ${Organ}_AT_4tunitIntersectNativeHMM_AlleleHMM.bed -b HMM_bed/${Organ}_PB6_HMM_*.bed \
| awk 'BEGIN {OFS="\t"}($4==$17){print $1,$2,$3,$4,$5,$6, $7, $8, $9, $10, $11, $12}' \
| awk 'BEGIN {OFS="\t"} (substr($4,1,1)=="P")  {print $0, "B6"} (substr($4,1,1)=="M")  {print $0, "CAST"}' >> ${Organ}_AT_4tunitIntersectNativeHMM_AlleleHMM_strain.bed

intersectBed -wb -s -a <(cat ${Organ}_AT_4tunitIntersectNativeHMM_AlleleHMM_strain.bed | cut -f 7-12) -b ${Organ}_AT_4tunitIntersectNativeHMM_AlleleHMM_strain.bed >${Organ}_AT_4tunitIntersectNativeHMM_intersectRegion_strain.bed