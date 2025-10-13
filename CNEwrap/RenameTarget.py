#!/usr/bin/env python3
def rename_maf(maffile,outmaf,refsp,tarsp,fmt="maf"):
	with open(maffile) as MAF, open(outmaf,'w') as OUTMAF:
		count_s=0
		for rec in MAF:
			if rec.startswith("s "):
				count_s+=1
				if not rec.startswith("s %s"%refsp) and not rec.startswith("s %s"%tarsp):
					if count_s%2==1:
						print(rec.replace("s ","s %s_"%refsp),file=OUTMAF,end="")
					if count_s%2==0:
						print(rec.replace("s ","s %s_"%tarsp),file=OUTMAF,end="")
				else:
					print(rec,file=OUTMAF,end="")
			else:
					print(rec,file=OUTMAF,end="")
if __name__=="__main__":
#information from bedfile
	maffile="target/ActMar/Tbai19.net.filter.axt.maf"
	outmaf="Tbai19.ActMar.net.filter.axt.maf"
	rename_maf(maffile,outmaf,"Tbai","ActMar",fmt="maf")
	#batch_write_align(bedfile,mafdir,aimsp,fmt="maf")
