#!/usr/bin/env python3
import os,sys
from glob import glob
from CNEwrap.utils import runcmd,runpool,get_system_id
from CNEwrap.TreeManipulation import traverse_from_leaf
class MZ(object):
	def __init__(self,treefile,refsp,suffix,AllMAF,cmdir):
		#command location
		ostype=get_system_id()
		if ostype=="ubuntu":
			subcmdir=os.path.join(cmdir,"x86_64")
		elif ostype=="centos":
			subcmdir=os.path.join(cmdir,"x86_64vc")
		else:
			raise ValueError("Error, can not determine your system!!")
			return(1)
		self.treefile=treefile
		self.refsp=refsp
		self.suffix=suffix
		self.AllMAF=AllMAF
		self.cmdir=subcmdir
		#self.treenodes = get_tree_node(self.treefile)
		self.treenodes = traverse_from_leaf(self.treefile,self.refsp)
		self.filename_map = {}
		self.species_number=len(self.treenodes)
		(self.finalmafs,self.commandlists)=self.batch_MergeMAF()
	def run_multiz(self,ABmaf,ACmaf):
		ABlist=ABmaf.replace("."+self.suffix,"").split('.')
		A=ABlist[0]
		B=".".join(ABlist[1:])
		inlevel=len(ABlist)
		if inlevel <= 2:
			self.filename_map[ABmaf]=ABmaf
		else:
			simname=A+".L"+str(inlevel)+"."+self.suffix
			assert self.filename_map[ABmaf] == simname, "%s and %s not equal"%(ABmaf,simname) 
		AClist=ACmaf.replace("."+self.suffix,"").split('.')
		A2=AClist[0]
		self.filename_map[ACmaf]=ACmaf
		assert A == A2, "the two input files are with different prefix"
		C=".".join(AClist[1:])
		ABCmaf=".".join([A,B,C,self.suffix])
		outlevel=len(ABlist)
		self.filename_map[ABCmaf]=A+".L"+str(inlevel+1)+"."+self.suffix
		command="{multiz} {ABmaf} {ACmaf} 1 > {ABCmaf}".format(
		multiz=os.path.join(self.cmdir,"multiz"),
		ABmaf=os.path.join(self.AllMAF,self.filename_map[ABmaf]),
		ACmaf=os.path.join(self.AllMAF,self.filename_map[ACmaf]),
		ABCmaf=os.path.join(self.AllMAF,self.filename_map[ABCmaf]))
		#whether to remove prervious produced intermediate maffile
		if inlevel > 2:
			rmmaf_cmd="rm {ABmaf}".format(ABmaf=os.path.join(self.AllMAF,self.filename_map[ABmaf]))
			command=command+" && "+ rmmaf_cmd
		return(ABCmaf,command)
	def MergeMAF(self,refmaf):
		maflist=[ ".".join([refmaf,qry,self.suffix]) for qry in self.treenodes ]
		commandlist=[]
		for i in range(len(maflist)-1):
			(maflist[i+1],command)=self.run_multiz(maflist[i],maflist[i+1])
			commandlist.append(command)
		return(maflist[i+1],commandlist)

	def batch_MergeMAF(self):
		self.treenodes.remove(self.refsp)
		AllMAF = os.path.abspath(self.AllMAF)
		mafs = glob("%s/*.maf"%self.AllMAF)
		refmafs = set([os.path.basename(maf).split('.')[0] for maf in mafs])
		finalmafs=[]
		commandlists=[]
		for refmaf in refmafs:
			(finalmaf,commandlist)=self.MergeMAF(refmaf)
			finalmafs.append(finalmaf)
			commandlists.append(commandlist)
		return(finalmafs,commandlists)
	def runList(self,commandlist):
		for command in commandlist:
			runcmd(command)
	def run_batch_MergeMAF(self):
		arglist=[(commandlist,) for commandlist in self.commandlists]
		runpool(self.runList,arglist,len(arglist))
	
	###do split and filter
	def FilterSplit(self,finalmaf,dryrun=False):
		os.makedirs("splitmaf",exist_ok=True)
		mafFilter=os.path.join(self.cmdir,"mafFilter")
		simfinalmaf=self.filename_map[finalmaf]
		filtermaf=simfinalmaf.replace(self.suffix,"FILTER."+self.suffix)
		Fcommand="{mafFilter} -minRow={species_number} -maxRow=500 -minCol=30 -minScore=20000 {simfinalmaf} > {filtermaf}".format(
			mafFilter=mafFilter,
			species_number=self.species_number,
			simfinalmaf=os.path.join(self.AllMAF,simfinalmaf),
			filtermaf=os.path.join(self.AllMAF,filtermaf))
		if dryrun:
			print(Fcommand)
		else:
			runcmd(Fcommand)
		mafSplit=os.path.join(self.cmdir,"mafSplit")
		Scommand="{mafSplit} -byTarget _.bed -useFullSequenceName splitmaf/{simfinalmaf} {filtermaf}".format(
			mafSplit=mafSplit,
			simfinalmaf=simfinalmaf.split(".")[0],
			filtermaf=os.path.join(self.AllMAF,filtermaf))
		if dryrun:
			print(Scommand)
		else:
			runcmd(Scommand)
	def batch_filter_split(self,ncpu=10,dryrun=False):
		argslist=[(finalmaf,dryrun) for finalmaf in self.finalmafs]
		maxcpu=ncpu if len(argslist) > ncpu else  len(argslist)
		runpool(self.FilterSplit,argslist,len(argslist))
	
if __name__ == "__main__":
	cmdir="/data/nfs/yancc/CNEWRAP/bin"
	suffix="net.filter.axt.maf"
	AllMAF="AllMAF"
	treefile="inputs/mammal.tre"
	refsp="simHuman"
	mz=MZ(treefile,refsp,suffix,AllMAF,cmdir)
	mz.run_batch_MergeMAF()
	mz.batch_filter_split()
	print(mz.filename_map)
