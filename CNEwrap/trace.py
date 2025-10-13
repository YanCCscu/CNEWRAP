#!/usr/bin/env python3
import pybedtools
from CNEwrap.MergeBed import parse_gerp_bed,parse_phast_bed,parse_gff,print_interval
from CNEwrap.ExtractMAFbyBed import batch_write_align
def trace_seqs(gerp_dir,phast_dir,splitmaf,gfffile,refgenome,bedkey="cne",fasdir="bed_fasta"):
	####MergeBed
	print('==='*10+"  parse gerp results   "+"==="*10)
	gerplist=parse_gerp_bed(gerp_dir)
	print('==='*10+"  parse phast results  "+"==="*10)
	phaslist=parse_phast_bed(phast_dir)
	gfflist=parse_gff(gfffile)
	print('==='*10+"    satrt mergebed     "+"==="*10)
	gerp_phas_bed=pybedtools.BedTool("\n".join(gerplist+phaslist),from_string=True)
	gffbed=pybedtools.BedTool("\n".join(gfflist),from_string=True)
	gerp_phas_nocds=gerp_phas_bed.sort().subtract(gffbed.sort())
	gerp_phas_nocds_rename=gerp_phas_nocds.sort().merge(c='4,5,6',o='collapse,distinct,distinct',delim="_")
	bedfile=print_interval(gerp_phas_nocds_rename,bedkey)
	print('==='*10+"   MergeBed succeed    "+"==="*10)
	
	####ExtractMAFbyBed
	batch_write_align(bedfile,splitmaf,refgenome,fasdir=fasdir)
	print("==="*10+"ExtractMAFbyBed succeed"+"==="*10)

if __name__ == "__main__":
	gerp_dir="gerp_dir"
	phast_dir="phast_dir"
	gfffile="inputs/simHuman.gff"
	
	bedfile="cne30.bed"
	refgenome="simHuman"
	fasdir="bed_fasta"
	splitmaf="splitmaf"
	trace_seqs(gerp_dir,phast_dir,splitmaf,gfffile,refgenome,bedfile=bedfile,fasdir="bed_fasta")
