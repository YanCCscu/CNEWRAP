#!/usr/bin/env python3
import sys
import os
import subprocess
import time
from glob import glob
from CNEwrap.utils import runcmd,runpool
from CNEwrap.TreeManipulation import get_tree_node

def get_prefix(sname):
	sname=os.path.basename(sname)
	return(".".join(sname.split('.')[:-1]))

def run_trf(fa,ref_dir,cmdir):
	trf_tmp_dir = os.path.join(refdir + "%s_trf.tmp"%(get_prefix(fa)))
	ref_trf=os.path.join(refdir + "ref_trf")
	os.makedirs(trf_tmp_dir, exist_ok=True)
	os.makedirs(ref_trf, exist_ok=True)
	command="{command} -bed -bedAt={ref_dir}/{ref_trf}/{out}.bed -trf={trfpath} -tempDir={trf_tmp_dir} {fa} {ref_dir}/{ref_trf}/{out}".format(
		command=os.path.join(cmdir,"trfBig"),
		ref_dir=ref_dir,
		out=get_prefix(fa),
		trf_tmp_dir=trf_tmp_dir,
		ref_trf=ref_trf,
		fa=fa)
	runcmd(command)
	command="{command} maskfasta -fi {fa} -bed {ref_dir}/{ref_trf}/{out}.bed -fo {fa} -soft".format(
		command=os.path.join(cmdir,"trfBig"),
		fa=fa,
		ref_dir=ref_dir,
		ref_trf=ref_trf,
		out=get_prefix(fa))
	rmcmd(command)	
	subprocess.run("/bin/rm -rf %s"%(trf_tmp_dir))

def fa2bit(fa,ref_dir,cmdir):
	command="{command} {fasta} {ref_dir}/{out}.2bit".format(command=os.path.join(cmdir,"faToTwoBit"),
		fasta=fa,
		ref_dir=ref_dir,
		out=get_prefix(fa))
	runcmd(command)

def fasize(fa,ref_dir,cmdir):
	command="{command} {fasta} -detailed > {ref_dir}/{out}.size".format(command=os.path.join(cmdir,"faSize"),
		fasta=fa,
		ref_dir=ref_dir,
		out=get_prefix(fa))
	runcmd(command)

def fasplit(ref_file,splitn,ref_dir,refgenome,cmdir):
	command="{fasplit} sequence {ref} {splitn} {refpath}".format(fasplit=os.path.join(cmdir,"faSplit"),
		ref=ref_file,
		splitn=splitn,
		refpath=os.path.join(ref_dir,refgenome))
	runcmd(command)

#genomedir: a dir containing all genome fasta files
#refgenome: set reference genome key abbreviation, like Pgut/Tele et al.,
#cmdir: path of bin will be used
def prepare_genomes(genomedir,refgenome,treefile,splitn,cmdir,dotrf=False):
	curdir = os.getcwd()
	ref_dir = os.path.join(curdir, "reference")
	tar_dir = os.path.join(curdir, "target")
	os.makedirs(ref_dir, exist_ok=True)
	os.makedirs(tar_dir, exist_ok=True)
	treenodes=get_tree_node(treefile)

	all_files =glob(os.path.join(genomedir, "*.fa*")) 
	ref_file = glob(os.path.join(genomedir, "%s*"%refgenome))
	tar_files = all_files[:]
	tar_files.remove(ref_file[0])
	suffix = tar_files[0].split('.')[-1]

	# Split fas into $plitn file and put them into directory ${refgenome}.cut
	print("prepare reference genome ...")
	fasplit(ref_file[0],splitn,ref_dir,refgenome,cmdir)
	#------------------
	refGP = glob(os.path.join(ref_dir, "%s*.fa*"%refgenome))
	print("\033[1;32msplit ref-genome into files: {0}\033[0m".format(refGP))
	splitn = len(refGP)
	GPlist=[(GP,ref_dir,cmdir) for GP in refGP]
	if dotrf == True:
		runpool(do_trf,GPlist,splitn)
	else:
		runpool(fasize,GPlist,splitn)
		runpool(fa2bit,GPlist,splitn)
	# Deal with query genomes
	print("prepare target genomes ...")
	tar_list=[]
	for g in tar_files:
		print("deal with %s ..."%g)
		tar_sp = os.path.basename(g).replace(f'.{suffix}', '')
		if tar_sp in treenodes:
			tar_dir_sp = os.path.join(tar_dir, tar_sp)
			os.makedirs(tar_dir_sp, exist_ok=True)
			tar_list.append((g,tar_dir_sp,cmdir))
	tarlen=len(tar_list)
	runpool(fasize,tar_list,tarlen)
	runpool(fa2bit,tar_list,tarlen)
	return(ref_dir,tar_dir)

if __name__ == '__main__':
	pydir=os.path.dirname(sys.argv[0])
	cmdir=os.path.join(pydir,"../bin")
	refgenome = "simHuman"
	genomedir = "mammals"
	treefile = "inputs/mammals.tre"
	refgenome=sys.argv[1]
	genomedir=sys.argv[2]
	treefile = sys.argv[3] 
	splitn = sys.argv[4]
	prepare_genomes(genomedir,refgenome,treefile,splitn,cmdir,False)
	#python CNEwrap/FormatGenome.py Tele input_genomes 20
