#!/usr/bin/env python3
import os,glob
from CNEwrap.MergeMAF import MZ
from CNEwrap.utils import runcmd,runpool
def remove_trailing_numbers(s):
    index = len(s)
    while index > 0 and s[index-1].isdigit():
        index -= 1
    return s[:index]

def rename_chrom(maffile,outmaf,refsp,tarsp):
        refsp=remove_trailing_numbers(refsp)
        with open(maffile) as MAF, open(outmaf,'w') as OUTMAF:
                count_s=0
                for rec in MAF:
                        if rec.startswith("s "):
                                count_s+=1
                                if not rec.startswith("s %s_"%refsp) and not rec.startswith("s %s_"%tarsp):
                                        if count_s%2==1:
                                                print(rec.replace("s ","s %s_"%refsp),file=OUTMAF,end="")
                                        if count_s%2==0:
                                                print(rec.replace("s ","s %s_"%tarsp),file=OUTMAF,end="")
                                else:
                                        print(rec,file=OUTMAF,end="")
                        else:
                                        print(rec,file=OUTMAF,end="")

def softlink_maf(workdir,suffix,targetdir="target",dryrun=False,increment=False):
	os.makedirs("AllMAF",exist_ok=True)
	#files = glob.glob('AllMAF/*')
	#for fl in files:
	#	os.remove(fl)
	for r,d,f in os.walk(os.path.join(workdir,targetdir),followlinks=True):
		for fname in f:
			if fname.endswith(suffix):
				refname=fname.split(".")[0]
				tagname=os.path.basename(r)
				source_name=os.path.join(r,fname)
				link_name=os.path.join(workdir,"AllMAF",refname+"."+tagname+"."+suffix)
				if dryrun:
					os.symlink(source_name,link_name)
					print("create softlink "+source_name+" --> "+link_name)
				else:
					if os.path.exists(link_name):
						if increment == False:
							os.remove(link_name)
							print("{link_name} exists,now removing it".format(link_name=link_name))
						else:
							continue
					os.symlink(source_name,link_name)
					print("create softlink "+source_name+" --> "+link_name)

def rename_maf(workdir,suffix,targetdir="target",dryrun=False,increment=False):
	os.makedirs("AllMAF",exist_ok=True)
	for r,d,f in os.walk(os.path.join(workdir,targetdir),followlinks=True):
		for fname in f:
			if fname.endswith(suffix):
				refname=fname.split(".")[0]
				tagname=os.path.basename(r)
				source_name=os.path.join(r,fname)
				link_name=os.path.join(workdir,"AllMAF",refname+"."+tagname+"."+suffix)
				if False:
					print("renamed chromosomes and move "+source_name+" --> "+link_name)
				else:
					if os.path.exists(link_name):
						if increment == False:
							os.remove(link_name)
							print("{link_name} exists,now removing it".format(link_name=link_name))
						else:
							continue
					rename_chrom(source_name,link_name,refname,tagname)
					print("renamed chromosomes and move "+source_name+" --> "+link_name)

def merge_maf(treefile,refgenome,suffix,cmdir,mafdir="AllMAF",targetdir="target",renamechr=False, \
	invoke_sge=False,sge_para="-l vf=10G,p=1",maxjobs=100,check_interval=10,bypass_rename=False,increment=False, dryrun=False):
	workdir=os.getcwd()
	if not bypass_rename:
		if renamechr:
			rename_maf(os.path.abspath(workdir),suffix,targetdir="target",dryrun=dryrun,increment=increment)
		else:		
			softlink_maf(os.path.abspath(workdir),suffix,targetdir="target",dryrun=dryrun,increment=increment)
	mz=MZ(treefile,refgenome,suffix,mafdir,cmdir)
	if invoke_sge:
		from CNEwrap.SubmitSGE import subjobs
		subjobs(mz.commandlists,maxjobs=maxjobs,check_interval=check_interval,jobname="mergemaf",sge_para=sge_para,dryrun=dryrun)
	else:
		if dryrun:
			for cmdline in mz.commandlists:
				print("\n".join(cmdline))
		else:
			mz.run_batch_MergeMAF()
	mz.batch_filter_split(dryrun=dryrun)

if __name__ == "__main__":
	cmdir="/home/ycc/OSshare/STAGEV/bin"
	suffix="net.filter.axt.maf"
	mafdir="AllMAF"
	treefile="inputs/mammal.tre"
	refgenome="simHuman"
	#softlink_maf(os.path.abspath(workdir),suffix)
	#mz=MZ(treefile,refsp,suffix,AllMAF,cmdir)
	#mz.batch_filter_split()
	#print(mz.finalmafs)
	merge_maf(treefile,refgenome,suffix,mafdir,cmdir)
