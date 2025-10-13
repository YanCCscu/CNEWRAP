#!/usr/bin/env python3
import os,sys
import numpy as np
from scipy.stats import gamma
from glob import glob
from multiprocessing import Process
from CNEwrap.GetPerID  import batch_write_perid
from CNEwrap.ForwardGenomics import make_pheno,merge_perid,make_anc_tree,run_FG 
from CNEwrap.EvoAcc import parse_dist,evo_acc
from CNEwrap.PhyloP import batch_phylop 
from CNEwrap.Phast import phylofit,cat_mafs 
from CNEwrap.RepMatNu import replace_matrix

def runFG(fasdir,mytree,aimsp,cmdir,outdir="FGout"):
	batch_write_perid(fasdir,mytree,cmdir)
	os.makedirs(outdir,exist_ok=True)
	(IDlist,peridlocal,peridglobal) = merge_perid(fasdir)
	pheno=make_pheno(mytree,aimsp)
	cne_anc_tre=make_anc_tree(mytree)
	run_FG(cne_anc_tre,pheno,IDlist,peridglobal,peridlocal,outdir,cmdir)

def cal_rep_mat(maffile,seqnumber=10000,minlen=500):
	#finalmaf_path="splitmaf/{finalmaf}".format(finalmaf=finalmaf.split(".")[0])
	#maffiles=glob("%s*.maf"%finalmaf_path)
	#for maffile in maffiles:
	print("Calculate dist matrix for "+maffile,file=sys.stderr)
	my_lom=replace_matrix(maffile,seqnumber,minlen)
	mafdist=maffile+'.dist'
	with open(mafdist,'w') as DistMX:                                                                                                                    
		print(my_lom.format(fmt="%6.3f"),file=DistMX) 
	return(mafdist)

def spevo(fasdir,species,treefile,bgfile=None,fraction=1,gaps=1,outdir="SPout",distfile=None):
	os.makedirs(outdir,exist_ok=True)
	SUMFILE=open("%s/SPsum.tsv"%outdir,'w')
	STATFILE=open("%s/statsum.txt"%outdir,"w")
	alidirlist=glob("%s/*"%fasdir)
	print("elementID\tsite\tJscore\tfgIC\tbgIC\tsiteIC\tdnadist\tICscore\tacc_score\tpvalue\tsignal\tbgAA\tfgAA",file=STATFILE)
	print("elementID\tseqlen\tvar_sites\tsig_sites\tsum_score\tmax_score\tave_score\tpvalue\tacc",file=SUMFILE)
	for alidir in alidirlist:
		evo_acc(alidir,species,treefile,distfile,bgfile=None,STATFILE=STATFILE,SUMFILE=SUMFILE,fraction=1,gaps=1)
	STATFILE.close()
	SUMFILE.close()

def phylop(fasdir,phylopmod,fgbranch,cmdir,mytree,phylop_dir="PPout",ncpu=20):
	alidirlist=glob("%s/*"%fasdir)
	os.makedirs(phylop_dir,exist_ok=True)
	phylout=os.path.join(phylop_dir,"PhyloP.txt")
	with open(phylout,'w') as PHYLOUT:
		header='elementID\tstart\tend\tname\tnull_scale\talt_scale\talt_subscale\tlnlratio\tpval'
		print(header,file=PHYLOUT)
		for alidir in alidirlist:
			if os.path.isdir(alidir) and len(os.listdir(alidir)) > 0:
				print("runing %s ... "%alidir)
				res=batch_phylop(alidir,phylopmod,fgbranch,cmdir,phylop_dir,ncpu)
				for r in sum(res,[]):
					print(r,end="",file=PHYLOUT)
def acc_cne(fasdir,treefile,fgspecies,phylopmod,cmdir,bgfile=None,fraction=1,gaps=1,evo_method="all",distfile=None):
	maflist=[maffile for maffile in glob(os.path.join("MAFBlock","*.maf"))]
	outmaf="mergedsplit.maf"
	if not phylopmod:
		print("Calculate phylofit based on MAFBlock")
		phylopmod=phylofit(maflist,treefile,cmdir,outmaf)
	if not distfile and evo_method in ["SP","SPPP","all"]:
		if not os.path.exists("mergedsplit.maf"):
			cat_mafs(maflist,outmaf)
		distfile=cal_rep_mat("mergedsplit.maf",1000,100)
	if evo_method == "FG":
		runFG(fasdir,treefile,fgspecies,cmdir,outdir="FGout")
	if evo_method == "SP":
		spevo(fasdir,fgspecies,treefile,bgfile=None,fraction=1,gaps=1,outdir="SPout",distfile=distfile)
	if evo_method == "PP":
		phylop(fasdir,phylopmod,fgspecies,cmdir,treefile,phylop_dir="PPout",ncpu=20)
	if evo_method == "SPPP":
		pw2=Process(target=spevo,args=(fasdir,fgspecies,treefile,None,1,1,"SPout",distfile))
		pw3=Process(target=phylop,args=(fasdir,phylopmod,fgspecies,cmdir,treefile))
		pw2.start()
		pw3.start()
		pw2.join()
		pw3.join()
	if evo_method == "all":
		pw1=Process(target=runFG,args=(fasdir,treefile,fgspecies,cmdir))
		pw2=Process(target=spevo,args=(fasdir,fgspecies,treefile,None,1,1,"SPout",distfile))
		pw3=Process(target=phylop,args=(fasdir,phylopmod,fgspecies,cmdir,treefile))
		pw1.start()
		pw2.start()
		pw3.start()
		pw1.join()
		pw2.join()
		pw3.join()

if __name__ == "__main__":
	fasdir="bed_fasta"
	treefile="CNE.tre"
	cmdir="/data/nfs/yancc/CNEWRAP/bin"
	fgspecies="Tbai"
	phylopmod="mergedsplit.mod"
	#acc_test(fasdir,mytree,aimsp,cmdir)
	#spspsite(fasdir,aimsp,bgfile=None,fraction=1,gaps=1)
	acc_cne(fasdir,treefile,fgspecies,cmdir,bgfile=None,fraction=1,gaps=1,evo_method="PP")

