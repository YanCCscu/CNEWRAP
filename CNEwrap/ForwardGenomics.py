#!/usr/bin/env python3
import sys,os
from glob import glob
from CNEwrap.TreeManipulation import mark_tree_node,get_tree_node 
from CNEwrap.utils import runcmd

g='\033[1;32m'
r='\033[1;91m'
b='\033[0m'

def run_FG(cne_anc_tree,pheno,IDlist,peridglobal,peridlocal,outdir,cmdir):
	FG=os.path.join(cmdir,"ForwardGenomics","forwardGenomics.R")
	print(FG)
	forward_command="{FG} \
--tree={cne_anc_tree} \
--elementIDs={IDlist} \
--expectedPerIDs={cmdir}/ForwardGenomics/lookUpData/expPercentID_CNE.txt \
--weights={cmdir}/ForwardGenomics/lookUpData/branchWeights_CNE.txt \
--method=all \
--minLosses=1 \
--listPheno={pheno} \
--globalPid={peridglobal} \
--localPid={peridlocal} \
--outFile={outdir}/Output.txt \
--outPath={outdir}/outpath \
--verbose=FALSE \
--thresholdConserved=0".format(FG=FG,
			cne_anc_tree=cne_anc_tree,
			IDlist=IDlist,
			cmdir=cmdir,
			pheno=pheno,
			peridglobal=peridglobal,
			peridlocal=peridlocal,
			outdir=outdir) 
	runcmd(forward_command)

def merge_perid(perid_dir):
	peridglobal_list=glob(os.path.join(perid_dir,"*/*.peridglobal"))
	peridlocal_list=glob(os.path.join(perid_dir,"*/*.peridlocal"))
	assert len(peridglobal_list) == len(peridlocal_list), "Number of peridglobal and number peridlocal is different"
	peridlocal="allcne.peridlocal"
	peridglobal="allcne.peridglobal"
	IDlist="ID.list"
	OUT_PL=open(peridlocal,'w')
	OUT_PG=open(peridglobal,'w')
	OUT_ID=open(IDlist,'w')
	global_title=""
	local_title=""
	for gf in peridglobal_list:
		with open(gf) as GF:
			if not global_title:
				global_title=GF.readline()
				print(global_title,end="",file=OUT_PG)
			else:
				GF.readline()
			for gl in GF:
				print(gl,end="",file=OUT_PG)
				gid=gl.split()[0]
				print(gid,file=OUT_ID)
	for lf in peridlocal_list:
		with open(lf) as LF:
			if not local_title:
				local_title=LF.readline()
				print(local_title,end="",file=OUT_PL)
			else:
				LF.readline()
			for ll in LF:
				print(ll,end="",file=OUT_PL)
	return(IDlist,peridlocal,peridglobal)

def make_pheno(CNEtre,aimsp):
	aimsplist=aimsp.strip().split(",")
	tree_nodes=get_tree_node(CNEtre)
	pheno=open("ID.pheno",'w')
	print("species pheno",file=pheno)
	for sp in tree_nodes:
		if sp in  aimsplist:
			print("%s 1"%sp,file=pheno)
		else:
			print("%s 0"%sp,file=pheno)
	return("ID.pheno")
def make_anc_tree(CNEtre):
	cne_anc_tre=mark_tree_node(CNEtre)
	treefile="CNEanc.tre"
	with open(treefile,'w') as TREE:
		print(cne_anc_tre,file=TREE,end="")
	return(treefile)
	
if __name__ == "__main__":
	cmdir="/data/nfs/yancc/CNE_TEST/STAGEII/bin"
	outdir="FGout"
	CNEtre="mammal.tre"
	aimsp="simHuman"
	perid_dir="bed_fasta"
	os.makedirs(outdir,exist_ok=True)
	(IDlist,peridlocal,peridglobal) = merge_perid(perid_dir)
	pheno=make_pheno(CNEtre,aimsp)
	cne_anc_tre=make_anc_tree(CNEtre)
	run_FG(cne_anc_tre,pheno,IDlist,peridglobal,peridlocal,outdir,cmdir)
