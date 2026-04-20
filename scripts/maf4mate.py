#!/usr/bin/env python3
import sys
import os
import subprocess
from collections import defaultdict
"""
Functions for  progress of extraction, convertion and concatenation
"""
def parse_names(line):
    """
    Parse NAMES=... line.
    Return contig_names and species_names.
    """
    names = line.split("=", 1)[1].strip().split(",")
    species = [n.split("_", 1)[0] for n in names]
    return names, species
"""
Extract 4D sites from PHAST SS file and concatenate by species (strict version).
"""
def ss2fasta(ss_file, out_fasta):
    # species -> list of bases
    species_seqs = defaultdict(list)
    contig_names = None
    species_names = None
    with open(ss_file) as f:
        for lineno, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if line.startswith("NAMES"):
                contig_names, species_names = parse_names(line)
                continue
            fields = line.split()
            if len(fields) < 5:
                continue
            comp3 = fields[3]
            if contig_names is None:
                raise RuntimeError("NAMES line not found before data")
            if len(comp3) != len(contig_names):
                raise ValueError(
                    f"Line {lineno}: comp3 length ({len(comp3)}) "
                    f"!= NAMES count ({len(contig_names)})"
                )
            per_species_bases = defaultdict(list)
            for base, species in zip(comp3, species_names):
                per_species_bases[species].append(base)
            for species, bases in per_species_bases.items():
                non_star = {b for b in bases if b != "*"}
                if len(non_star) == 0:
                   non_star={"-"} 
                if len(non_star) > 1:
                    raise ValueError(
                        f"Line {lineno}: species {species} has "
                        f"multiple bases at one 4D site: {sorted(non_star)}"
                    )
                species_seqs[species].append(non_star.pop())

    with open(out_fasta, "w") as out:
        for species in sorted(species_seqs):
            seq = "".join(species_seqs[species])
            out.write(f">{species}\n")
            for i in range(0, len(seq), 60):
                out.write(seq[i:i+60] + "\n")

"""
Extract 4D sites from a MAF file using PHAST msa_view.
"""
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
    print(f" ----- produce gffout {output_gff} for {scaffold}")
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
def maf24d(maf,gff="Tbai.cds.gff",outdir="maf4d"):
    msa_view = os.path.join(
        script_dir(),
        "..", "bin", "phast-1.3", "bin", "msa_view"
    )
    if not os.path.exists(outdir): 
        os.makedirs(outdir, exist_ok=True)
    if not os.path.isfile(msa_view):
        raise ValueError(f"{msa_view} not found!")
    # 1. Extract scaffold
    scaffold = extract_scaffold(maf)
    print(f"Processing scaffold: {scaffold}")
    # 2. Filter GFF
    gff_out = os.path.join(outdir, f"{scaffold}.gff")
    print(f"filter_gff: {gff} {gff_out} {scaffold}")
    filter_gff(gff, scaffold, gff_out)
    # 3. Extract 4D sites (SS format)
    ss_out = os.path.join(outdir, f"{scaffold}.4d.ss")
    run_msa_view(maf,gff_out,ss_out,msa_view)
    print(f"produce gffout {ss_out}")
    # 4. Convert SS → FASTA
    fas_out = os.path.join(outdir, f"{scaffold}.4d.fas")
    ss2fasta(ss_out,fas_out)
    print(f"produce gffout {fas_out}")
    return fas_out

"""
Merge multiple FASTA files by concatenating sequences with the same title.
"""
def read_fasta(path):
    """
    Read FASTA file.
    Return dict: title -> sequence
    """
    seqs = {}
    title = None
    buf = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if title:
                    seqs[title] = "".join(buf)
                title = line[1:].strip().split("_")[0]
                buf = []
            else:
                buf.append(line)
        if title:
            seqs[title] = "".join(buf)
    return seqs

def read_maf(path):
    """
    read MAF (Multiple Alignment Format) file
    """
    seqs = defaultdict(list)
    with open(path) as f:
        for line in f:
            if line.startswith("s "):  # MAF line starts with 's ' 
                parts = line.strip().split()
                if len(parts) < 7:
                    continue
                src = parts[1]
                title = src.split("_")[0]
                sequence = parts[6]
                seqs[title].append(sequence)
    
    return {title: "".join(fragments) for title, fragments in seqs.items()}

def concate_4d(out_fasta, input_files):
    merged = defaultdict(list)
    for filepath in input_files:
        if not os.path.exists(filepath):
            print(f"{filepath} not exists, skip ...", file=sys.stderr)
            continue
        ext = os.path.splitext(filepath)[1].lower()
        if ext in ['.fas', '.fa', '.fasta']:
            seqs = read_fasta(filepath)
        elif ext == '.maf':
            seqs = read_maf(filepath)
        else:
            print(f"{ext} of {filepath} can not tell, use fasta or maf file please", file=sys.stderr)
            seqs = read_fasta(filepath)
        for title, seq in seqs.items():
            merged[title].append(seq)
    #write outfile
    with open(out_fasta, "w") as out:
        for title in sorted(merged):
            full_seq = "".join(merged[title])
            out.write(f">{title}\n")
            for i in range(0, len(full_seq), 60):
                out.write(full_seq[i:i+60] + "\n")

'''
if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(f"Usage: {sys.argv[0]} input.ss output.fasta")
    ss2fasta(sys.argv[1], sys.argv[2])
    args = parse_args()
    maf = args.maf
    gff = args.gff
    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)
    maf24d(maf,gff,outdir)
    if len(sys.argv) < 3:
        sys.exit(
            f"Usage: {sys.argv[0]} output.fasta input1.fas input2.fas ..."
        )
    concate_4d(sys.argv[1], sys.argv[2:])
'''
