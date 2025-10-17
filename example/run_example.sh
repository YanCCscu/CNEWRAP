#!/bin/bash
. /home/ycc/anaconda3/bin/activate /home/ycc/cnewrap_env 
CNEHOME=$(realpath ../)
export PATH=$PATH:${CNEHOME}/bin

#make random tree if you do not have a species tree
#python ../scripts/mk_random_tree.py sp1 sp2 sp3 ...

#produce lastz.dist
python  $CNEHOME/CNEwrap/TreeManipulation.py CNE.tre Tbai

gunzip  input_genomes/*.gz
#you can modify lastz.dist to give suitable distance of each species pair [far,medium,near]
python $CNEHOME/cnewrap.py align -r Tbai -t CNE.tre -g input_genomes -n 20 -d lastz.dist 
# ro you can run the above scripts in SGE clusters with --invoke_sge

python $CNEHOME/cnewrap.py merge -r Tbai -t CNE.tre -m AllMAF -s net.filter.axt.maf 
# ro you can run the above scripts in SGE clusters with --invoke_sge

python $CNEHOME/cnewrap.py scne -r Tbai -t CNE.tre -b MAFBlock -p phast_dir -m splitmaf -g gerp_dir 
# ro you can run the above scripts in SGE clusters with --invoke_sge

python $CNEHOME/cnewrap.py trace -r Tbai -a Tbai.cds.gff -e tbaicne -p phast_dir -s splitmaf -f bed_fasta -g gerp_dir 

python $CNEHOME/cnewrap.py evolve -f Tbai -t CNE.tre -d bed_fasta -p mergedsplit.mod --evo_method all

#refer ../scripts/run_phyloacc_single.sh for phyloacc method
#rm -rf allcne.peridglobal allcne.peridlocal AllMAF bed_fasta CNEanc.tre FGout gerp_dir ID.list ID.pheno lastz.dist MAFBlock mafinfo.table mergedsplit.log mergedsplit.maf mergedsplit.maf.dist mergedsplit.mod phast_dir PPout reference splitmaf SPout target tbaicne30.bed 
