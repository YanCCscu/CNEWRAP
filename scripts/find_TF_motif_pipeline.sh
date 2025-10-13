#Download TF.motif.txt from CisBP
#http://cisbp.ccbr.utoronto.ca/index.php
########################################################################################################
#Homo Sapiens' files in /data/nfs/tcy/AP_2020/find_motif/homo_sapiens

#use PoSSuM for mapping TFs on genome fasta
perl /data/nfs/OriginTools/find_motif/PoSSuM-1_3/bin/Transfac2PoSSuM-PSSM.pl $1 > PoSSuM-PSSM.motifs.txt     #$1 is motifs' file path,like /data/nfs/tcy/AP_2020/find_motif/homo_sapiens/pwms_all_motifs/; $2 is genome.fa file
/data/nfs/OriginTools/find_motif/PoSSuM-1_3/bin/possumfreqs -db $2 -dna > ${2}.motif.frequencies.txt
/data/nfs/OriginTools/find_motif/PoSSuM-1_3/bin/possumdist -pr PoSSuM-PSSM.motifs.txt -dna -freq ${2}.motif.frequencies.txt -pdis dist.gz
/data/nfs/OriginTools/find_motif/PoSSuM-1_3/bin/possumsearch -db $2 -pr PoSSuM-PSSM.motifs.txt  -dna -pval 1e-05 -simple  -pdis dist.gz 1> ${2}.TF.motifs.txt
perl /data/nfs/OriginTools/find_motif/PoSSuM-1_3/bin/get_PoSSuM_results.pl ${2}.TF.motifs.txt > ${2}.TF.motifs.bed

#if focus on specific TFs
awk '{if($7~"FOXD3") print}' ../find_motif/homo_sapiens/TF_Information_all_motifs_plus.txt > TF_Information_all_motifs_plus_FOXD3.txt     #grep FOXD3's information
awk '{print $4}' TF_Information_all_motifs_plus_FOXD3.txt > FOXD3_human_motif_ID     #make FOXD3 motif ID list
for i in `cat FOXD3_human_motif_ID`;do awk '{if($8~"'$i'") print}' ${2}.TF.motifs.bed;done >> FOXD3_human_motif_mapping
