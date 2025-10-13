#!/usr/bin/env python3
import sys,os
from glob import glob
from ete3 import Tree
from Bio import SeqIO
from CNEwrap.utils import runcmd,runpool
import subprocess
g='\033[1;32m'
r='\033[1;91m'
b='\033[0m'
MINLEN=30

def GetPerID(Seq,Seqanc):
	assert len(Seq)==len(Seqanc),\
	"seq and seq_anc not in same length, make sure they were aligned ...\n"
	total_len=0
	matchbp=0
	for a,b in zip(Seq,Seqanc):
		if a == b:
			if a == '-':
				continue
			else:
				matchbp+=1
		total_len+=1
	return (matchbp,total_len)
		
#parse tree return three relative dicts		
def ParseTree(treefile):
	mytree=open(treefile)
	while True:
			#skip empty lines
			treestr=mytree.readline()
			if treestr != "":
				break
	t=Tree(treestr,format=1)
	nodesAnc=[]
	leafRoot=[]
	ancnode=""
	ancnum2name={}
	for node in t.traverse("preorder"):
		if not node.is_root():
			node.name=node.name.split("_")[0]
			nodesAnc.append((node.name,node.get_ancestors()[0].name))
		else:
			ancnode=node.name
	for node in t.traverse("preorder"):
		if node.is_leaf():
			leafRoot.append((node.name,ancnode))
			

	for node in t.traverse("postorder"):
		if not node.is_leaf():
			ancnum=node.name
			node.name=node.get_children()[0].name.split("-")[0]+\
			"-"+node.get_children()[1].name.split("-")[0]
			ancnum2name[ancnum]=node.name
	return (leafRoot,nodesAnc,ancnum2name)

def run_prank(myfas,mytree,cmdir):
	prank=os.path.join(cmdir,"prank")
	myprefix=myfas.replace('.fas','').replace('.fa','')
	prankcommand="%s -d=%s -o=%s -t=%s -showtree -showanc -keep -prunetree -prunedata -quiet -seed=10"%(prank,myfas,myprefix,mytree)
	myancfas=myprefix+".anc.fas"
	myanctree=myprefix+".anc.dnd"
	try:
		runcmd(prankcommand,runstdout=subprocess.DEVNULL, runstderr=subprocess.DEVNULL)
		#runcmd(prankcommand)
		return(myancfas,myanctree)
	except:
		return(myancfas,myanctree)
def write_perid(myfas,mytree,cmdir):
	myprefix=myfas.replace('.fas','').replace('.fa','')
	myancfas=myprefix+".anc.fas"
	myanctree=myprefix+".anc.dnd"
	SeqID=os.path.basename(myfas).replace('.fas','').replace('.fa','')
	outdir=os.path.dirname(myfas) if myfas.find('/')!=-1 else os.getcwd()
	if (not os.path.exists(myancfas)) or (not os.path.exists(myanctree)):
		(myancfas,myanctree)=run_prank(myfas,mytree,cmdir)
		if (not os.path.exists(myancfas)) or (not os.path.exists(myanctree)):
			print("PASS: %s not produced to will be ignored ..."%(myancfas))
			return(2)
	fasdict={ rec.id.split("_")[0] : rec.seq for rec in SeqIO.parse(myancfas,'fasta') }
	(leafRoot,nodesAnc,ancnum2name)=ParseTree(myanctree)
	#------------parse tree and make dict for
	# leaves : ancestor for peridglobal
	# internode :ancestor for peridlocal
	# internode : intername for peridlocal id

#------------write peridglobal
	print(g,"%s output leaves identity refer to root node in dir: %s ..."%(SeqID,outdir),b)
	with open(outdir+"/"+SeqID+".peridglobal",'w') as PERIDGLO:
		global_ids=[];leaves_sps=[]
		for (m,n) in sorted(leafRoot, key=lambda x: x[0]):
			(match,alllen)=GetPerID(fasdict[m],fasdict[n])
			leaves_sps.append(m)
			if alllen>=MINLEN:
				percent_id=match/alllen
				global_ids.append(str(percent_id))
			else:
				global_ids.append("NA")
		print("species"," "," ".join(leaves_sps),file=PERIDGLO,sep="")
		print(SeqID," "," ".join(global_ids),file=PERIDGLO,sep="")
#-------------write peridlocal
	print(g,"%s output internode identity refer to its MRCA in dir: %s ..."%(SeqID,outdir),b)
	(leafRoot,nodesAnc,ancnum2name)=ParseTree(myanctree)
	with open(outdir+"/"+SeqID+".peridlocal",'w') as PERIDLOC:
		print("branch id pid",file=PERIDLOC)
		for (k,p) in nodesAnc:
			(match,alllen)=GetPerID(fasdict[k],fasdict[p])
			if alllen>=MINLEN:
				if k in ancnum2name:
					print("%s %s %f"%(ancnum2name[k],SeqID,match/alllen),file=PERIDLOC,sep="") 
				else:
					print("%s %s %f"%(k,SeqID,match/alllen),file=PERIDLOC,sep="")
			else:
				print("%s %s NA"%(k,SeqID),file=PERIDLOC,sep="")
def batch_write_perid(fasdir,mytree,cmdir,parallel=False):
		faslist=glob(os.path.join(fasdir,"*/*.fas"))
		argslist=[(myfas,mytree,cmdir) for myfas in faslist if "anc.fas" not in myfas]
		if parallel:
			runpool(write_perid,argslist,len(argslist))
		else:
			for myfas in faslist:
				if "anc.fas" not in myfas:
					write_perid(myfas,mytree,cmdir)
		


###################################################################################################
####---------main process------####################################################################
###################################################################################################
if __name__ == "__main__":
	myfas="bed_fasta/Tbai11Tbai_hic_scaffold_12.maf/CNE00026060.fas"
	fasdir="bed_fasta"
	mytree="CNE.tre"
	cmdir="/data/nfs/yancc/CNEWRAP/bin"
	#write_perid(myfas,mytree,cmdir)
	batch_write_perid(fasdir,mytree,cmdir)
