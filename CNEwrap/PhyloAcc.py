#!/usr/bin/env python3
import os
from glob import glob
from Bio import AlignIO,SeqIO
from CNEwrap.utils import runcmd,runcmd_pipe,runpool
def runphyloacc(bedfas,phylopmod,fgbranch,PhyloP,phylop_dir):
	fasbase=os.path.basename(bedfas).replace(".fas","")
	phylout=os.path.join(phylop_dir,fasbase+".phylop")
	fasbed=os.path.join(phylop_dir,fasbase+".bed")
	with open(fasbed,'w') as BED:
		rec=next(SeqIO.parse(bedfas,"fasta"))
		pchr=rec.id
		pend=len(rec.seq)
		BED.write("{chr}\t{start}\t{end}\n".format(chr=pchr,start=0,end=pend))
	phyloacc_command='phyloacc.py -a {bedfas} -b {fasbed} -m {phylopmod} -t "{fgbranch}" --local && rm {fasbed}'.format(
		fgbranch=fgbranch.replace(",",";"),
		fasbed=fasbed,
		phylopmod=phylopmod,
		bedfas=bedfas
	)
	#with open(phylout,'w') as PHYLOUT:
	outstd=runcmd_pipe(phylop_command)
	return(outstd.readlines()[1:])
	
def batch_phyloacc(fasdir,phylopmod,fgbranch,cmdir,phylop_dir="PPout",ncpu=20):
	os.makedirs(phylop_dir,exist_ok=True)
	phyloP=os.path.join(cmdir,"phast-1.3/bin/phyloP")
	alilist=glob("%s/*.fas"%fasdir)
	alilist=[ bedfas for bedfas in alilist if not ".anc.fas" in bedfas]
	argslist=[(bedfas,phylopmod,fgbranch,phyloP,phylop_dir) for bedfas in alilist]
	maxcpu=ncpu if len(argslist) > ncpu else len(argslist) 
	res=runpool(runphylop,argslist,maxcpu)
	return(res)
if __name__ == "__main__":
	cmdir="/data/nfs/yancc/CNEWRAP/bin"
	phyloFit=os.path.join(cmdir,"phast-1.3/bin/phyloFit")
	phyloP=os.path.join(cmdir,"phast-1.3/bin/phyloP")
	fasdir="bed_fasta/Tbai19Tbai_hic_scaffold_2840.maf"
	phylopmod="mergedsplit.mod"
	fgbranch='Tbai'
	phylop_dir="PPout"
	os.makedirs(phylop_dir,exist_ok=True)
	header='#chr\tstart\tend\tname\tnull_scale\talt_scale\talt_subscale\tlnlratio\tpval'
	print(header)
	res=batch_phylop(fasdir,phylopmod,fgbranch,cmdir,phylop_dir="PPout",ncpu=20)
	for r in sum(res,[]):
		print(r,end="")
	#cmdlist=batch_phylop_sge(fasdir,phylopmod,fgbranch,cmdir,phylop_dir="PPout")
	#print(cmdlist)
