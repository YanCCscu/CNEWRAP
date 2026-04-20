#!/usr/bin/env python3
"""
Extract 4D sites from PHAST SS file and concatenate by species (strict version).
"""
import sys
from collections import defaultdict
def parse_names(line):
    """
    Parse NAMES=... line.
    Return contig_names and species_names.
    """
    names = line.split("=", 1)[1].strip().split(",")
    species = [n.split("_", 1)[0] for n in names]
    return names, species
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
if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(f"Usage: {sys.argv[0]} input.ss output.fasta")
    ss2fasta(sys.argv[1], sys.argv[2])
