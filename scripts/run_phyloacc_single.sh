#!/bin/bash
source /data/share/OriginTools/OriginTools/anaconda3/bin/activate phyloacc
CNEHOME=/data/nfs/yancc/CNEWRAP
find  bed_fasta -name "*.fas"|grep -v anc.fas >fasta.list
grep \> $(head -1 fasta.list)|sed 's/>//' >samples.list
python $CNEHOME/scripts/concatenate_bed.py -a samples.list -f fasta.list -o concatenated
$CNEHOME/bin/tree_doctor --name-ancestors mergedsplit.mod >mergelabel.mod
casespA="Rkio"
case=ACC
phyloacc.py -a concatenated.fas \
-b concatenated.bed \
-m mergelabel.mod \
-t "$casespA" -o phyloacc_st_${case}M0 --local --overwrite
`grep snakemake phyloacc_st_${case}M0/phyloacc_st_${case}M0.log |rev | cut -d' ' -f2- | rev`
phyloacc_post.py -i phyloacc_st_${case}M0 -o phyloacc_st_${case}M0 --overwrite
