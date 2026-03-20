# obtain the annotation information of CNEs for each of the two reference genomes
grep "^Tbai" ../Tbai_inputs/bed_fasta/*maf/*.gff |sed 's/.\+://' >Tbai_refer_Tbai_cnes.gff
grep "^Tbai" ../Tele_inputs/bed_fasta/*maf/*.gff |sed 's/.\+://' >Tele_refer_Tbai_cnes.gff

# fix and unify into the Tbai coordinate format
awk 'BEGIN{OFS="\t"} {if($2 > $3){temp=$2; $2=$3; $3=temp} print $0}' Tbai_refer_Tbai_cnes.gff > Tbai_refer_Tbai_fixed.gff
awk 'BEGIN{OFS="\t"} {if($2 > $3){temp=$2; $2=$3; $3=temp} print $0}' Tele_refer_Tbai_cnes.gff > Tele_refer_Tbai_fixed.gff

# obtain the common CNEs
bedtools intersect -a Tbai_refer_Tbai_fixed.gff -b Tele_refer_Tbai_fixed.gff -g Tbai_genome.size -wo >CompareCNEs_Tbai.table

# obtain the uniq CNEs
bedtools subtract -a Tbai_refer_Tbai_fixed.gff -b Tele_refer_Tbai_fixed.gff -A >Tbai-Tele_Tbai_diffr.gff
bedtools subtract -a Tele_refer_Tbai_fixed.gff -b Tbai_refer_Tbai_fixed.gff -A >Tele-Tbai_Tbai_diffr.gff

# merge all the common and uniq CNEs
cat CompareCNEs_Tbai.table Tbai-Tele_Tbai_diffr.gff Tele-Tbai_Tbai_diffr.gff >merged_TbaiCNEs.gff
