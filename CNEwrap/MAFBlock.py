#!/usr/bin/env python3
import os
from Bio import AlignIO
from glob import glob
from CNEwrap.utils import runcmd,runpool
def splitperblock(maffile,blockdir):
	count=0
	mafinfo=[]
	mafbase=os.path.basename(maffile).replace(".maf","")
	blockmafdir=os.path.join(blockdir,mafbase)
	os.makedirs(blockmafdir, exist_ok=True)
	for alignment in AlignIO.parse(maffile, "maf"):
		scaf=""
		start=0
		refchr=alignment[0].id
		outmaf=os.path.join(blockmafdir,"%s.%07d"%(refchr,count)+".maf")
		for rec in alignment:
			if refchr in rec.id:
				scaf=rec.id
				start=int(rec.annotations['start'])
			if "_" in rec.id:
				rec.id=rec.id.split("_")[0]
			else:
				rec.id=rec.id.split(".")[0]
		#print("\rDealing with: %s"%outmaf,end="")
		AlignIO.write(alignment,outmaf,"maf")
		count+=1
		mafinfo.append((outmaf,scaf,start))
	return(mafinfo)

def batch_splitblock(mafdir,blockdir):
	maflist=glob(os.path.join(mafdir,"*.maf"))
	argslist=[(maffile,blockdir) for maffile in maflist]
	nthread= len(argslist) if len(argslist) <= 40 else 40
	res=runpool(splitperblock,argslist,nthread)
	return(res)


if __name__ == "__main__":
	from pprint import pprint
	os.makedirs("MAFBlock",exist_ok=True)
	res=batch_splitblock("splitmaf","MAFBlock")
	pprint(res)
