#!/usr/bin/env python3
import os
from glob import glob
from Bio import AlignIO
from CNEwrap.utils import runcmd,runpool

def cat_mafs(inmaflist,outmaf="mergedsplit.maf"):
	with open(outmaf, 'w', encoding='utf-8') as f_out:
		for maffile in inmaflist:
			with open(maffile, 'r', encoding='utf-8') as f_in:
				f_out.write(f_in.read())
	return(outmaf)

def phylofit(maflist,treefile,cmdir,outmaf="mergedsplit.maf"):
	phyloFit=os.path.join(cmdir,"phast-1.3/bin/phyloFit")
	phyloFitpara="--subst-mod REV --EM --precision HIGH --msa-format MAF"
	cat_mafs(maflist,outmaf)
	mafbase=os.path.basename(outmaf).replace(".maf","")
	mafmod=os.path.basename(outmaf).replace(".maf",".mod")
	maflog=os.path.basename(outmaf).replace(".maf",".log")
	phylofit_command="{phyloFit} {para} --tree {treefile} --out-root {mafout} -l {maflog} {maffile}".format(
		para=phyloFitpara,
		phyloFit=phyloFit,
		treefile=treefile,
		mafout=mafbase, #output will auto add .mod
		maflog=maflog,
		maffile=outmaf)
	runcmd(phylofit_command)
	return(mafmod)

def runphast_sge(outmaf,phastmod,scaf,start,phast_dir,phastCons):
	commandlist=[]
	last_dir = os.path.basename(os.path.dirname(outmaf))
	os.makedirs(os.path.join(phast_dir,last_dir),exist_ok=True)
	mafbase=os.path.basename(outmaf).replace(".maf","")
	mafbed=os.path.join(phast_dir,last_dir,mafbase+".con.bed")
	phastlog=os.path.join(phast_dir,last_dir,mafbase+".phastCons.log")
	phastwig=os.path.join(phast_dir,last_dir,mafbase+".score.wig")
	phastConspara="--expected-length 45 --target-coverage 0.3 --rho 0.3"
	phastcons_command="{phastCons} {para} --most-conserved {mafbed} --log {phastlog} {maffile} {phastmod} > {phastwig}".format(
		phastCons=phastCons,
		para=phastConspara,
		mafbed=mafbed,
		phastlog=phastlog,
		maffile=outmaf,
		phastmod=phastmod,
		phastwig=phastwig) 
	rm_tmp='/bin/rm %s %s'%(phastwig,phastlog)
	rename_chr="[[ ! -s %s ]] && /bin/rm %s || sed -i 's/^\S\+/%s/' %s"%(mafbed,mafbed,scaf,mafbed)
	commandlist=[phastcons_command,rm_tmp,rename_chr]
	return commandlist

#SGE mode only return commands for sge submition
def batch_phast_sge(mafinfo,phastmod,phast_dir,cmdir):
	os.makedirs(phast_dir,exist_ok=True)
	#phyloFit=os.path.join(cmdir,"phast-1.3/bin/phyloFit")
	phastCons=os.path.join(cmdir,"phast-1.3/bin/phastCons")
	commandlist=[runphast_sge(outmaf,phastmod,scaf,start,phast_dir,phastCons) for (outmaf,scaf,start) in mafinfo]
	return commandlist

def runphast(outmaf,phastmod,scaf,start,phast_dir,phastCons):
	last_dir = os.path.basename(os.path.dirname(outmaf))
	os.makedirs(os.path.join(phast_dir,last_dir),exist_ok=True)
	mafbase=os.path.basename(outmaf).replace(".maf","")
	mafout=os.path.join(phast_dir,last_dir,mafbase+".phast")
	mafbed=os.path.join(phast_dir,last_dir,mafbase+".con.bed")
	phastlog=os.path.join(phast_dir,last_dir,mafbase+".phastCons.log")
	phastwig=os.path.join(phast_dir,last_dir,mafbase+".score.wig")
	phastConspara="--expected-length 45 --target-coverage 0.3 --rho 0.3"
	phastcons_command="{phastCons} {para} --most-conserved {mafbed} --log {phastlog} {maffile} {phastmod} > {phastwig}".format(
		phastCons=phastCons,
		para=phastConspara,
		mafbed=mafbed,
		phastlog=phastlog,
		maffile=outmaf,
		phastmod=phastmod,
		phastwig=phastwig) 
	runcmd(phastcons_command)
	#rename_chr="sed -i 's/^\S\+/%s/' %s"%(scaf,mafbed)
	rename_chr="[[ ! -s %s ]] && /bin/rm %s || sed -i 's/^\S\+/%s/' %s"%(mafbed,mafbed,scaf,mafbed)
	runcmd(rename_chr)
	os.remove(phastwig)
	os.remove(phastlog)
	
def batch_phast(mafinfo,phastmod,phast_dir,cmdir,ncpu=40):
	os.makedirs(phast_dir,exist_ok=True)
	phastCons=os.path.join(cmdir,"phast-1.3/bin/phastCons")
	argslist=[(outmaf,phastmod,scaf,start,phast_dir,phastCons) for (outmaf,scaf,start) in mafinfo]
	maxcpu=ncpu if len(argslist) > ncpu else  len(argslist) 
	runpool(runphast,argslist,maxcpu)

if __name__ == "__main__":
	mafdir="MAFBlock_test"
	treefile="CNE.tre"
	outdir="phast_dir_test"
	cmdir="/data/nfs/yancc/CNEWRAP/bin"
	refsp='Tbai'
	phyloFit=os.path.join(cmdir,"phast-1.3/bin/phyloFit")
	phastCons=os.path.join(cmdir,"phast-1.3/bin/phastCons")
	alignmentlist=[]
	mafinfo=[]
	maflist=[maffile for maffile in glob(os.path.join(mafdir,"*.maf"))]
	mafmod=phylofit(maflist,treefile,cmdir)
	mafinfo=[(outmaf,"Tbai",1) for outmaf in maflist] 	
	#runphast(maffile,treefile,outdir,phyloFit,phastCons)
	batch_phast(mafinfo,mafmod,outdir,cmdir)
