#!/usr/bin/env python3
import os
from CNEwrap.FormatGenome import prepare_genomes 
from CNEwrap.RunLastz import LZ
from CNEwrap.utils import get_system_id,print_nested_list,runcmd,runpool


def batch_lz_dry(referenceDir,query2bit,divtime,cmdir,lzstep=1):
	mylz=LZ(referenceDir,query2bit,divtime,cmdir,lzstep)
	commandlist=mylz.batch_lastz_chain2maf()
	return commandlist

def runList(commandlist):
	for command in commandlist:
		runcmd(command)

def genome_align(genomedir,refgenome,treefile,splitn,divtime,cmdir,invoke_sge=False,sge_para="-l vf=40G,p=4 -pe smp 4",maxjobs=100,check_interval=10,dryrun=False):
	ostype=get_system_id()
	if ostype=="ubuntu":
		subcmdir=os.path.join(cmdir,"x86_64")
	elif ostype=="centos":
		subcmdir=os.path.join(cmdir,"x86_64vc")
	else:
		raise ValueError("Error, can not determine your system!!")
		return(1)
	(ref_dir,tar_dir)=prepare_genomes(genomedir,refgenome,treefile,splitn,subcmdir,False)
	tar_sps=os.listdir(tar_dir)
	divdict={}
	if not divtime in ["far","medium","near"]:
		assert os.path.exists(divtime),"make sure the lastz.dist file exist"
		with open(divtime) as DIV:
			for dt in DIV:
				dtl=dt.strip().split()
				assert dtl[0] == refgenome, "make sure the reference genome %s not exists"%dtl[0]
				assert dtl[1] in tar_sps, "%s not exists in target species"%dtl[1]
				divdict[dtl[1]]=dtl[2]
	else:
		for tarsp in tar_sps:
			divdict[tarsp]=divtime
			
	if invoke_sge:
		from CNEwrap.SubmitSGE import subjobs
		scripts=[]
		for tarsp in tar_sps:
			qury2bit=os.path.join(tar_dir,tarsp,tarsp+".2bit")
			commandlists=batch_lz_dry(ref_dir,qury2bit,divdict[tarsp],cmdir)
			for commandlist in commandlists:
				scripts.append(commandlist)
		subjobs(scripts,maxjobs=maxjobs,check_interval=check_interval,jobname="lzalign",sge_para=sge_para,dryrun=dryrun)
		
	else:
		scripts=[]
		for tarsp in tar_sps:
			qury2bit=os.path.join(tar_dir,tarsp,tarsp+".2bit")
			if dryrun:
				commandlists=batch_lz_dry(ref_dir,qury2bit,divdict[tarsp],cmdir)
				print_nested_list(commandlists)
			else:
				commandlists=batch_lz_dry(ref_dir,qury2bit,divdict[tarsp],cmdir)
				for commandlist in commandlists:
					scripts.append((commandlist,))
		ncpu=len(scripts) if len(scripts)<maxjobs else maxjobs
		runpool(runList,scripts,ncpu)

if __name__ == '__main__':
	refgenome = "simHuman"
	genomedir = "inputs/mammals"
	treefile="inputs/mammal.tre"
	cmdir="../bin"
	splitn = 2
	divtime="near"
	genome_align(genomedir,refgenome,treefile,splitn,divtime,cmdir,invoke_sge=False,sge_para="-l vf=4G,p=2 -pe smp 2",maxjobs=100,check_interval=10,dryrun=False)
