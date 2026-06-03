#!/usr/bin/env python3
import sys
import os
import subprocess
import time
from glob import glob
from CNEwrap.utils import runcmd,runpool

"""
COMMENT:
align every target sequences to reference with lastz
EXAMPLE:python CNEwrap/RunLastz.py reference target/simCow/simCow.2bit /data/nfs/yancc/CNE_TEST/STAGEII/bin
"""
class LZ(object):
	def __init__(self,referenceDir,query2bit,divtime,cmdir):
		self.divtime=divtime
		self.referenceDir=referenceDir
		self.query2bit=query2bit
		self.cmdir=cmdir
		self.para=self.get_para()
	def get_para(self):
		lastz_dist=os.path.join(self.cmdir,"lastz_dist")
		if self.divtime == "near":
			linearGap="{lastz_dist}/param/near.linearGap".format(lastz_dist=lastz_dist)
			chainPar="-minScore=5000 -linearGap={linearGap}".format(linearGap=linearGap)
			lastzPar="E=150 H=2000 K=4500 L=2200 M=254 O=600 Q={lastz_dist}/human-chimp.v2.q T=2 Y=15000".format(lastz_dist=lastz_dist)
		elif self.divtime == "medium":
			linearGap="medium"
			chainPar="-minScore=3000 -linearGap={linearGap}".format(linearGap=linearGap)
			lastzPar="E=30 H=2000 K=3000 L=2200 M=50 O=400 Q={lastz_dist}/HoxD70.q T=1 Y=9400".format(lastz_dist=lastz_dist)
		elif self.divtime == "far":
			linearGap="loose"
			chainPar="-minScore=5000 -linearGap={linearGap}".format(linearGap=linearGap)
			lastzPar="E=30 H=2000 K=2200 L=6000 M=50 O=400 Q={lastz_dist}/HoxD55.q T=2 Y=3400".format(lastz_dist=lastz_dist)
		else:
			print("Unrecognized a divergence time parameter!!")
			sys.exit(3)
		return(linearGap,chainPar,lastzPar)
	def lastz(self,ref2bit,lastzPar,axtout):
		lastz_D=os.path.join(self.cmdir,"lastz-distrib-1.03.54/bin/lastz_D")
		command="{lastz_D} {ref2bit}[multiple] {query2bit} --chain --gapped --notransition --format=axt {lastzPar} > {axtout}".format(
				lastz_D=lastz_D,
				ref2bit=ref2bit,
				query2bit=self.query2bit,
				lastzPar=lastzPar,
				axtout=axtout)
		runcmd(command)
				
	def batch_lastz(self):
		print("#"*20+"\nrun align target to reference with lastz in batch...\n"+"#"*20)
		querydir=os.path.dirname(self.query2bit)
		querybase=os.path.basename(self.query2bit.replace(".2bit",""))
		ref2bits=glob(os.path.join(self.referenceDir, "*.2bit"))
		lastzPar=self.get_para()
		runargs=[]
		for ref2bit in ref2bits:
			refbase=os.path.basename(ref2bit)
			outfile=os.path.join(querydir,querybase+"_"+refbase+".axt")
			runargs.append((ref2bit,lastzPar[2],outfile))
		runpool(self.lastz,runargs,len(runargs))

	def chain2maf(self,axtfile,runPar,ref2bit):
		#command location
		axtChain = os.path.join(self.cmdir,"axtChain")
		chainSort = os.path.join(self.cmdir,"chainSort")
		chainNet = os.path.join(self.cmdir,"GenomeAlignmentTools/bin/chainNet")
		NetFilter = os.path.join(self.cmdir,"GenomeAlignmentTools/bin/NetFilterNonNested.perl")
		netToAxt = os.path.join(self.cmdir,"netToAxt")
		axtToMaf = os.path.join(self.cmdir,"axtToMaf")
		#file rername
		ref_path_base=ref2bit.replace(".2bit","")
		ref_base=os.path.basename(ref_path_base)
		qry_path_base=self.query2bit.replace(".2bit","")
		qry_base=os.path.basename(qry_path_base)
	
		chain_outbase=axtfile.replace(".axt","")
		qrydir=os.path.dirname(axtfile)
	
		fref_size=ref_path_base+".size"
		fqry_size=qry_path_base+".size"
		faxtchain=chain_outbase+".axt.chain"
		fsortchain=chain_outbase+".axt.chain.sorted"
		frefnet=os.path.join(qrydir,ref_base+".net")
		fqrynet=os.path.join(qrydir,qry_base+"_"+ref_base+".net")
		frefnetF=frefnet+".filter"
		fqrynetF=fqrynet+".filter"
		fnetaxt=frefnetF+".axt"
		faxtmaf=fnetaxt+".maf"

		axt2chain="{axtChain} {chainPar} {axtfile} {refsp} {qrysp} {axtchain}".format(
			axtChain=axtChain,
			chainPar=runPar[1],
			axtfile=axtfile,
			refsp=ref2bit,
			qrysp=self.query2bit,
			axtchain=faxtchain)
	
		sortchain="{chainSort} {axtchain} {sortchain}".format(
			chainSort=chainSort,
			axtchain=faxtchain,
			sortchain=fsortchain)

		chainnet="{chainNet} -linearGap={linearGap} -rescore -tNibDir={refsp} -qNibDir={qrysp} {sortchain} {refsp_size} {qrysp_size} {refnet} {qrynet}".format(
			chainNet=chainNet,
			linearGap=runPar[0],
			refsp=ref2bit,
			qrysp=self.query2bit,
			sortchain=fsortchain,
			refsp_size=fref_size,
			qrysp_size=fqry_size,
			refnet=frefnet,
			qrynet=fqrynet)

		netfilter="{NetFilter} -minScore 20000 -minSizeT 4000 -minSizeQ 4000 {refnet} > {refnetF}".format(
			NetFilter=NetFilter,
			refnet=frefnet,
			refnetF=frefnetF)

		net2axt="{netToAxt} {refnetF} {axtchain} {refsp} {qrysp} {netaxt}".format(
			netToAxt=netToAxt,
			refnetF=frefnetF,
			axtchain=faxtchain,
			refsp=ref2bit,
			qrysp=self.query2bit,
			netaxt=fnetaxt)

		axt2maf="{axtToMaf} {netaxt} {refsp_size} {qrysp_size} {axtmaf}".format(
			axtToMaf=axtToMaf,
			netaxt=fnetaxt,
			refsp_size=fref_size,
			qrysp_size=fqry_size,
			axtmaf=faxtmaf)
		runcmd(axt2chain)
		runcmd(sortchain)
		runcmd(chainnet)
		runcmd(netfilter)
		runcmd(net2axt)
		runcmd(axt2maf)
	
	def batch_chain2maf(self):
		print("#"*30+"\nconvert and filter axt to maf in batch...\n"+"#"*30)
		querydir=os.path.dirname(self.query2bit)
		querybase=os.path.basename(self.query2bit.replace(".2bit",""))
		ref2bits=glob(os.path.join(self.referenceDir, "*.2bit"))
		runPar=self.get_para()
		runargs=[]
		#(self,axtfile,chainPar,ref2bit)
		for ref2bit in ref2bits:
			refbase=os.path.basename(ref2bit)
			axtfile=os.path.join(querydir,querybase+"_"+refbase+".axt")
			runargs.append((axtfile,runPar,ref2bit))
		runpool(self.chain2maf,runargs,len(runargs))

if __name__ == '__main__':
	import argparse
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument('-r', action='store', dest="referenceDir", help='reference genome for aligning to',required=True)
	parser.add_argument('-q', action='store', dest="query2bit", help='Directory containing all genome sequence files')
	align_choices=["far","medium","near"]
	parser.add_argument('-d', action='store', dest="divtime", choices=align_choices,default='medium',help='distance between species candidate: near medium far')
	parser.add_argument('-c', action='store', dest="cmdir", help='path of command directory')
	args = parser.parse_args()
	#referenceDir = "reference"
	#query2bit = "target/simCow/simCow.2bit"
	#divtime="near"
	#cmdir="/data/nfs/yancc/CNE_TEST/STAGEII/bin"
	
	mylz=LZ(args.referenceDir,args.query2bit,args.divtime,args.cmdir)
	mylz.batch_lastz()
	mylz.batch_chain2maf()
