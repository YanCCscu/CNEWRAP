#!/usr/bin/env python3
import sys
import os
import subprocess
import time
from glob import glob
from CNEwrap.utils import runcmd,runpool,get_system_id

"""
COMMENT:
align every target sequences to reference with lastz
EXAMPLE:python CNEwrap/RunLastz.py reference target/simCow/simCow.2bit /data/nfs/yancc/CNE_TEST/STAGEII/bin
"""
class LZ(object):
	def __init__(self,referenceDir,query2bit,divtime,cmdir,lzstep=1):
		self.divtime=divtime
		self.referenceDir=referenceDir
		self.query2bit=query2bit
		self.cmdir=cmdir
		self.lzstep=lzstep
		self.para=self.get_para()
	def get_para(self):
		lastz_dist=os.path.join(self.cmdir,"lastz_dist")
		if self.divtime == "near":
			linearGap="{lastz_dist}/param/near.linearGap".format(lastz_dist=lastz_dist)
			chainPar="-minScore=5000 -linearGap={linearGap}".format(linearGap=linearGap)
			lastzPar="--step={step} E=150 H=2000 K=4500 L=2200 M=254 O=600 Q={lastz_dist}/human-chimp.v2.q T=2 Y=15000".format(
				step=self.lzstep,lastz_dist=lastz_dist)
		elif self.divtime == "medium":
			linearGap="medium"
			chainPar="-minScore=3000 -linearGap={linearGap}".format(linearGap=linearGap)
			lastzPar="--step={step} E=30 H=2000 K=3000 L=2200 M=50 O=400 Q={lastz_dist}/HoxD70.q T=1 Y=9400".format(
				step=self.lzstep,lastz_dist=lastz_dist)
		elif self.divtime == "far":
			linearGap="loose"
			chainPar="-minScore=5000 -linearGap={linearGap}".format(linearGap=linearGap)
			lastzPar="--step={step} E=30 H=2000 K=2200 L=6000 M=50 O=400 Q={lastz_dist}/HoxD55.q T=2 Y=3400".format(
				step=self.lzstep,lastz_dist=lastz_dist)
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
		return command
				

	def chain2maf(self,axtfile,chainPar,linearGap,ref2bit):
		#command location
		ostype=get_system_id()
		if ostype=="ubuntu":
			subcmdir=os.path.join(self.cmdir,"x86_64")
		elif ostype=="centos":
			subcmdir=os.path.join(self.cmdir,"x86_64vc")
		else:
			raise ValueError("Error, can not determine your system!!")
			return(1)
		axtChain = os.path.join(self.cmdir,subcmdir,"axtChain")
		chainSort = os.path.join(self.cmdir,subcmdir,"chainSort")
		chainNet = os.path.join(self.cmdir,subcmdir,"chainNet")
		NetFilter = os.path.join(self.cmdir,subcmdir,"NetFilterNonNested.perl")
		netToAxt = os.path.join(self.cmdir,subcmdir,"netToAxt")
		axtToMaf = os.path.join(self.cmdir,subcmdir,"axtToMaf")
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
			chainPar=chainPar,
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
			linearGap=linearGap,
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
		return [axt2chain,sortchain,chainnet,netfilter,net2axt,axt2maf]
	
	def batch_lastz_chain2maf(self):
		querydir=os.path.dirname(self.query2bit)
		querybase=os.path.basename(self.query2bit.replace(".2bit",""))
		ref2bits=glob(os.path.join(self.referenceDir, "*.2bit"))
		(linearGap,chainPar,lastzPar)=self.get_para()
		commandlist=[]
		for ref2bit in ref2bits:
			refbase=os.path.basename(ref2bit)
			axtfile=os.path.join(querydir,querybase+"_"+refbase+".axt")
			lastz_command=self.lastz(ref2bit,lastzPar,axtfile)
			chain2maf_command=self.chain2maf(axtfile,chainPar,linearGap,ref2bit)
			command=[lastz_command]+chain2maf_command
			commandlist.append(command)
		return commandlist
	def runList(self,commandlist):
		for command in commandlist:
			runcmd(command)
	def run_batch_lastz_chain2maf(self):
		print("#"*20+"\nrun align target to reference with lastz in batch...\n"+"#"*20)
		argslist=[(commandlist,) for commandlist in self.batch_lastz_chain2maf()]
		runpool(self.runList,argslist,len(argslist))

if __name__ == '__main__':
	import argparse
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument('-r', action='store', dest="referenceDir", help='reference genome for aligning to',required=True)
	parser.add_argument('-q', action='store', dest="query2bit", help='Directory containing all genome sequence files')
	align_choices=["far","medium","near"]
	parser.add_argument('-d', action='store', dest="divtime", choices=align_choices,default='medium',help='distance between species candidate: near medium far')
	parser.add_argument('-s', action='store', dest="lzstep", type=int,default=1,help='step parameter for lastz')
	parser.add_argument('-c', action='store', dest="cmdir", help='path of command directory')
	args = parser.parse_args()
	#referenceDir = "reference"
	#query2bit = "target/simCow/simCow.2bit"
	#divtime="near"
	#cmdir="/data/nfs/yancc/CNE_TEST/STAGEII/bin"
	
	mylz=LZ(args.referenceDir,args.query2bit,args.divtime,args.cmdir,args.lzstep)
	mylz.run_batch_lastz_chain2maf()
