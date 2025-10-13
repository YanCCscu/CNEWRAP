#!/usr/bin/env python3
import os
from glob import glob
from CNEwrap.utils import runcmd,runpool
from Bio import AlignIO

def run_gerp(outmaf,scaf,start,treefile,refsp,cmdir,gerp_dir="gerp_dir"):
	gerpcol=os.path.join(cmdir,"gerpcol")
	gerpelem=os.path.join(cmdir,"gerpelem")

	command_gerpcol="{gerpcol} -t {treefile} -f {outmaf} -x .GERP.rates -e {refsp}".format(
	gerpcol=gerpcol,
	treefile=treefile,
	outmaf=outmaf,
	refsp=refsp)
	
	command_gerpelem="{gerpelem} -c {scaf} -s {start} -x .elems -w .ex -f {outmaf}.GERP.rates -d 0.01".format(
	gerpelem=gerpelem,
	scaf=scaf,
	start=start,
	outmaf=outmaf)
	runcmd(command_gerpcol)
	runcmd(command_gerpelem)
	try:
		os.remove(outmaf+".GERP.rates")
		os.remove(outmaf+".GERP.rates.ex")
	except:
		pass
	last_dir = os.path.basename(os.path.dirname(outmaf))
	os.makedirs(os.path.join(gerp_dir,last_dir),exist_ok=True)
	os.rename(outmaf+".GERP.rates.elems","{gerp_dir}/{last_dir}/{outmaf}.GERP.rates.elems".format(
	gerp_dir=gerp_dir,
	last_dir =last_dir, 
	outmaf=os.path.basename(outmaf)))

def batch_gerp(mafinfo,treefile,refsp,cmdir,gerp_dir="gerp_dir",ncpu=40):
	os.makedirs(gerp_dir,exist_ok=True)
	arglist=[(outmaf,scaf,start,treefile,refsp,cmdir,gerp_dir) for (outmaf,scaf,start) in mafinfo ]
	maxcpu= ncpu if len(arglist) > ncpu else len(arglist)
	runpool(run_gerp,arglist,maxcpu)

def run_gerp_sge(outmaf,scaf,start,treefile,refsp,cmdir,gerp_dir="gerp_dir"):
	gerpcol=os.path.join(cmdir,"gerpcol")
	gerpelem=os.path.join(cmdir,"gerpelem")
	last_dir = os.path.basename(os.path.dirname(outmaf))
	os.makedirs(os.path.join(gerp_dir,last_dir),exist_ok=True)
	command_gerpcol="{gerpcol} -t {treefile} -f {outmaf} -x .GERP.rates -e {refsp}".format(
	gerpcol=gerpcol,
	treefile=treefile,
	outmaf=outmaf,
	refsp=refsp)
	
	command_gerpelem="{gerpelem} -c {scaf} -s {start} -x .elems -w .ex -f {outmaf}.GERP.rates -d 0.01".format(
	gerpelem=gerpelem,
	scaf=scaf,
	start=start,
	outmaf=outmaf)
	rm_tmp='/bin/rm {outmaf}.GERP.rates;/bin/rm {outmaf}.GERP.rates.ex; mv {outmaf}.GERP.rates.elems {gerp_dir}/{last_dir}'.format(
	outmaf=outmaf,
	gerp_dir=gerp_dir,
	last_dir=last_dir)
	return [command_gerpcol,command_gerpelem,rm_tmp]

def batch_gerp_sge(mafinfo,treefile,refsp,cmdir,gerp_dir="gerp_dir"):
	os.makedirs(gerp_dir,exist_ok=True) 
	commandlist=[run_gerp_sge(outmaf,scaf,start,treefile,refsp,cmdir,gerp_dir) for (outmaf,scaf,start) in mafinfo ]
	return(commandlist)
#awk -v id=$blockid -v no=$blockNO -v OFS="\t" '{\$7=id".maf."no".GERP";printf("%s\t%s\t%s\t%s\t0\t+\n",\$1,\$2,\$3,\$7)}' $curdir/${mafblock}.c.GERP.rates.elems > $curdir/${mafblock}.c.GERP.rates.elems.f 

if __name__ == "__main__":
	maffile="MAFBlock/simHuman1chr7.simHuman.chr7.0000.maf.F"
	mafdir="MAFBlock_test"
	treefile="CNE.tre"
	refsp="Tbai"
	cmdir="/data/nfs/yancc/CNEWRAP/bin"
	mafinfo=[]
	for maffile in glob(os.path.join(mafdir,"*.maf")):
		mafinfo.append((maffile,'Tbai',1))
	#batch_gerp(mafinfo,treefile,refsp,cmdir)
	cmdlist=batch_gerp_sge(mafinfo,treefile,refsp,cmdir)
	print(cmdlist)
