#!/usr/bin/env python3
import os
from glob import glob
from multiprocessing import Process
from CNEwrap.MAFBlock import batch_splitblock
from CNEwrap.GERP import batch_gerp,batch_gerp_sge
from CNEwrap.Phast import phylofit,batch_phast,batch_phast_sge 


def cne_scan(mafblock,treefile,refgenome,splitmaf,phastdir,gerpdir,cmdir,mafinfofile="",phylopmod="",invoke_sge=False,ncpu=40,maxjobs=100,check_interval=1,sge_para="-l vf=1G,p=1",cne_method="both"):
	os.makedirs(mafblock,exist_ok=True)
	if not mafinfofile:
		mafinfo=batch_splitblock(splitmaf,mafblock)
		mafinfo=sum(mafinfo,[])
		with open("mafinfo.table","w") as MAFINFO:
			for mi in mafinfo:
				#print(mi,file=MAFINFO)
				print("%s\t%s\t%s"%mi,end="\n",file=MAFINFO)
		maflist=[maffile for (maffile,scaf,start) in mafinfo]
	else:
		maflist=[]
		mafinfo=[]
		with open(mafinfofile) as MAFFILE:
			for rec in MAFFILE:
				reclist=rec.strip().split()
				mafinfo.append((reclist[0],reclist[1],int(reclist[2])))
				maflist.append(reclist[0])
	if not phylopmod:
		mafmod=phylofit(maflist,treefile,cmdir)
	else:
		mafmod=phylopmod
	#for maffile in glob(os.path.join(mafblock,"*.maf")):
	#	(outmaf,scaf,start)=reshape_maf(maffile,refgenome)
	#	mafinfo.append((outmaf,scaf,start))
	if invoke_sge:
		from CNEwrap.SubmitSGE import subjobs
		gerp_commandlist=batch_gerp_sge(mafinfo,treefile,refgenome,cmdir,gerpdir)
		pw1=Process(target=subjobs,args=(gerp_commandlist,maxjobs,check_interval,"rungerp",sge_para))
		
		phast_commandlist=batch_phast_sge(mafinfo,mafmod,phastdir,cmdir)
		pw2=Process(target=subjobs,args=(phast_commandlist,maxjobs,check_interval,"runphast",sge_para))
		
	else:
		pw1=Process(target=batch_gerp,args=(mafinfo,treefile,refgenome,cmdir,gerpdir,ncpu))
	
		pw2=Process(target=batch_phast,args=(mafinfo,mafmod,phastdir,cmdir,ncpu))
	if cne_method == "GERP":
		pw1.start()
		pw1.join()
	if cne_method == "Phast":
		pw2.start()
		pw2.join()
	if cne_method == "both":
		pw1.start()
		pw2.start()
		pw1.join()
		pw2.join()

if __name__ == "__main__":
	splitmaf="splitmaf_test"
	mafblock="MAFBlock_test"
	treefile="CNE.tre"
	refgenome="Tbai"
	phastdir="phast_dir_test"
	cmdir="/data/nfs/yancc/CNEWRAP/bin"
	cne_scan(mafblock,treefile,refgenome,splitmaf,phastdir,cmdir,invoke_sge=True,maxjobs=10,check_interval=1)
