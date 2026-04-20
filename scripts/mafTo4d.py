#!/bin/env python3
import os
import sys
import glob
import argparse
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from maf4mate import maf24d,concate_4d
def parse_args():
	parser = argparse.ArgumentParser(
		description="Extract 4D sites from a MAF alignment"
	)
	parser.add_argument("maf", help="Input MAF file or directory")
	parser.add_argument("gff", help="Annotation GFF file")
	parser.add_argument(
		"-o", "--outdir",
		default="maf4d",
		help="Output directory (default: maf4d)"
	)
	parser.add_argument(
		"-c", "--cat4d",
		default="cat_4d.fas",
		help="Output concatenated fasta file (default: cat_4d.fas)"
	)
	return parser.parse_args()
args = parse_args()
mafs = args.maf
gff_file = args.gff
outdir = args.outdir
catfile = args.cat4d
if os.path.isdir(mafs):
	maf_files = glob.glob(f"{mafs}/*.maf")
else:
	maf_files=mafs
if not maf_files:
	print(f"Can not find any .maf files in {mafs}")
	sys.exit(1)
print(f"Start to convert {len(maf_files)} files...")

with ProcessPoolExecutor(max_workers=10) as executor:
	worker = partial(maf24d, gff=gff_file, outdir=outdir)
	executor.map(worker, maf_files)
processed_files = glob.glob(f"{outdir}/*.4d.fas")
if processed_files:
	print("concatenate 4d.fas files...")
	concate_4d(catfile,processed_files)
	print(f"Done! output saved to : {catfile}")
	print(f"The following command can be used to constructed phylogenetic tree with {catfile}")
	print("--"*10)
	print(f"trimal -in {catfile} -out {catfile.replace('.fas','')}.trim.fas -gt 0.9")
	print(f"iqtree -s {catfile.replace('.fas','')}.trim.fas -nt 20 -st DNA -bb 1000 -alrt 1000")
	print("--"*10)
else:
	print("Can not find any .4d.fas files,skip")
