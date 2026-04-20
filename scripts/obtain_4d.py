#!/usr/bin/env python3
"""
Extract 4D sites from a MAF file using PHAST msa_view.
"""
from ss2fas import ss2fasta
import argparse
import os
import subprocess
import sys
def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract 4D sites from a MAF alignment"
    )
    parser.add_argument("maf", help="Input MAF file")
    parser.add_argument("gff", help="Annotation GFF file")
    parser.add_argument(
        "-o", "--outdir",
        default="maf4d",
        help="Output directory (default: maf4d)"
    )
    return parser.parse_args()
def script_dir():
    """Return directory where this script is located."""
    return os.path.dirname(os.path.abspath(__file__))
def extract_scaffold(maf):
    """
    Extract scaffold name from the first 's' line in a MAF file.
    """
    with open(maf) as fh:
        for line in fh:
            if line.startswith("s"):
                return line.split()[1]

    raise RuntimeError(f"No scaffold found in MAF file: {maf}")
def filter_gff(gff, scaffold, output_gff):
    """
    Keep only GFF entries for the given scaffold.
    """
    with open(gff) as fin, open(output_gff, "w") as fout:
        for line in fin:
            if line.startswith(scaffold + "\t"):
                fout.write(line)
def run_msa_view(maf,gff,ss_out,msa_view):
    """
    Run msa_view and redirect stdout to ss_out.
    """
    args=[maf,"--4d","--features", gff]
    cmd = [msa_view] + args
    with open(ss_out, "w") as fout:
        subprocess.run(
            cmd,
            stdout=fout,
            stderr=sys.stderr,
            check=True
        )
def maf24d(maf,gff,outdir):
    msa_view = os.path.join(
        script_dir(),
        "..", "bin", "phast-1.3", "bin", "msa_view"
    )
    if not os.path.isfile(msa_view):
        raise ValueError(f"{msa_view} not found!")
    # 1. Extract scaffold
    scaffold = extract_scaffold(maf)
    print(f"Processing scaffold: {scaffold}")
    # 2. Filter GFF
    gff_out = os.path.join(outdir, f"{scaffold}.gff")
    filter_gff(gff, scaffold, gff_out)
    # 3. Extract 4D sites (SS format)
    ss_out = os.path.join(outdir, f"{scaffold}.4d.ss")
    run_msa_view(maf,gff_out,ss_out,msa_view)
    # 4. Convert SS → FASTA
    fas_out = os.path.join(outdir, f"{scaffold}.4d.fas")
    ss2fasta(ss_out,fas_out)
    return fas_out
if __name__ == "__main__":
    args = parse_args()
    maf = args.maf
    gff = args.gff
    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)
    maf24d(maf,gff,outdir)
