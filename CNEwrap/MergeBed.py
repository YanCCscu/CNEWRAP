#!/usr/bin/env python3
import os
#from glob import glob
from pathlib import Path
import pybedtools
MIN_LEN=30
def parse_gerp_bed(gerp_dir):
	#gerp_elems=glob(os.path.join(gerp_dir,"**","*.GERP.rates.elems"),recursive=True)
	GD=Path(gerp_dir)
	gerp_elems=list(GD.rglob("*.GERP.rates.elems"))
	bedstrings=[]
	N=1
	for elems in gerp_elems:
		with open(elems) as gbed:
			for line in gbed:
				llist=line.strip().split("\t")
				chrom=llist[0]
				start=llist[1]
				end=llist[2]
				cneid="GERP_%s_%04d"%(chrom,N)
				N+=1
				bedline="{chrom}\t{start}\t{end}\t{cneid}\t0\t+".format(
				chrom=chrom,
				start=start,
				end=end,
				cneid=cneid)
				bedstrings.append(bedline)
	return(bedstrings)
def parse_phast_bed(phast_dir):
	#phast_beds=glob(os.path.join(phast_dir,"*con.bed"))
	PD=Path(phast_dir)
	phast_beds=list(PD.rglob("*con.bed"))
	bedstrings=[]
	J=1
	for bed in phast_beds:
		with open(bed) as pbed:
			for line in pbed:
				llist=line.strip().split("\t")
				chrom=llist[0]
				#chrom=".".join(chromosome.split(".")[1:])
				start=llist[1]
				end=llist[2]
				cneid="Phas_%s_%04d"%(chrom,J)
				J+=1
				bedline="{chrom}\t{start}\t{end}\t{cneid}\t0\t+".format(
				chrom=chrom,
				start=start,
				end=end,
				cneid=cneid)
				bedstrings.append(bedline)
	return(bedstrings)
def parse_gff(gff_file):
	bedstrings=[]
	K=1
	with open(gff_file) as GFF:
		for f in GFF:
			if f.startswith("#"):
				pass
			else:
				fl=f.strip().split("\t")
				ftype=fl[2]
				if ftype == "CDS" or ftype == "cds":
					chrom=fl[0]
					start=int(fl[3])-1
					end=fl[4]
					cdsid="CDS%06d"%(K)
					K+=1
					strand=fl[7]
					bedline="{chrom}\t{start}\t{end}\t{cdsid}\t0\t+".format(
					chrom=chrom,
					start=start,
					end=end,
					cdsid=cdsid)
					bedstrings.append(bedline)
		return(bedstrings)
def print_interval(mybed,keyword="cne"):
	bedfile="{key}{wlen}.bed".format(key=keyword,wlen=MIN_LEN)
	BED=open(bedfile,'w')
	uniqid=[]
	cneNO=0
	for bed in mybed:
		cneNO+=1
		if bed.length <= MIN_LEN:
			continue
		elif bed.name in uniqid:
			uniqid.append(bed.name)
			bed.name=bed.name+".%d"%uniqid.count(bed.name)
		else:
			uniqid.append(bed.name)
		pbed=bed.fields
		cnename="CNE%08d"%cneNO
		print("\t".join(pbed[0:3]+[cnename]+pbed[4:]+[pbed[3]]),end="\n",file=BED)
	BED.close()
	return(bedfile)

if __name__ == "__main__":
	gerp_dir="MAFBlock"
	phast_dir="phast_dir_test"
	gff_file="simHuman.gff"
	gerplist=parse_gerp_bed(gerp_dir)
	phaslist=parse_phast_bed(phast_dir)
	gfflist=parse_gff(gff_file)
	gerp_phas_bed=pybedtools.BedTool("\n".join(gerplist+phaslist),from_string=True)
	gffbed=pybedtools.BedTool("\n".join(gfflist),from_string=True)
	#print(gerp_phas_bed.sort().intersect(gffbed.sort()))
	gerp_phas_nocds=gerp_phas_bed.sort().subtract(gffbed.sort())
	gerp_phas_nocds_rename=gerp_phas_nocds.sort().merge(c='4,5,6',o='collapse,distinct,distinct',delim="_")
	print_interval(gerp_phas_nocds_rename)

