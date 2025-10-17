<img src="https://private-user-images.githubusercontent.com/42508147/502349357-e5f6e716-10be-4563-b7ec-087e2ecdde57.png?jwt=eyJ0eXAiOiJKV1QiLCJhbGciOiJIUzI1NiJ9.eyJpc3MiOiJnaXRodWIuY29tIiwiYXVkIjoicmF3LmdpdGh1YnVzZXJjb250ZW50LmNvbSIsImtleSI6ImtleTUiLCJleHAiOjE3NjA2NjYzODMsIm5iZiI6MTc2MDY2NjA4MywicGF0aCI6Ii80MjUwODE0Ny81MDIzNDkzNTctZTVmNmU3MTYtMTBiZS00NTYzLWI3ZWMtMDg3ZTJlY2RkZTU3LnBuZz9YLUFtei1BbGdvcml0aG09QVdTNC1ITUFDLVNIQTI1NiZYLUFtei1DcmVkZW50aWFsPUFLSUFWQ09EWUxTQTUzUFFLNFpBJTJGMjAyNTEwMTclMkZ1cy1lYXN0LTElMkZzMyUyRmF3czRfcmVxdWVzdCZYLUFtei1EYXRlPTIwMjUxMDE3VDAxNTQ0M1omWC1BbXotRXhwaXJlcz0zMDAmWC1BbXotU2lnbmF0dXJlPWYxNWQ5MjkzZTFmNDg0ZmFjMGY5NDRlMDYwNzU3OWYwNTA5NDhjOTQyZmYxNGMxNGE1NmExNzdmYTg5NWZjZjMmWC1BbXotU2lnbmVkSGVhZGVycz1ob3N0In0.FPEYWxK77lpYcOyM8R-LT4fdLHkYBAlvKRIhrg38udI" alt="LOGO of CNEwrap" width="4" height="4" />
# CNEwrap: A scalable toolkit with a novel algorithm for large-scale genome-wide accelerated conserved non-coding elements detection
---
**CNEwrap** is a Python-based pipeline designed to identify **conserved noncoding elements (CNEs)** across multiple genomes, perform multiple sequence alignments, and conduct evolutionary analyses, specifically focusing on CNEs.

## Install
```bash
git clone https://github.com/YanCCscu/CNEWRAP.git
conda env create -f cnewrap_env.yml
conda activate cnewrap_env
````

Alternatively, you can install with the Conda command:

```bash
conda install yccscucib::cnewrap
```

-----

## Usage

### Input files

#### Mandatory:

1.  **Genome Directory:** A directory **containing** genome sequence files in FASTA format. The prefix of each file name will be used as the **species identifier** in the provided species tree.
2.  **Species Tree:** A phylogenetic tree in **Newick format** that includes branch length information. (Note: If you do not have a species tree, you can construct one using a provided script—see `obtain_4d.sh`—after the alignment process).
3.  **Reference GFF:** A GFF file for the reference genome. **Only CDS items** within this file will be used for masking or analysis.

#### Optional:

1.  **`lastz.dist`:** Contains distance information between the reference species and every non-reference species.
2.  **`.mod`:** A phylogenetic model file used for the **phast** process.
3.  **`.dist`:** A distance file used for the **EvoAcc** process.

### Basic Command

```bash
python $CNEHOME/cnewrap.py -h
````
---

CNEwrap provides a series of modular subcommands that can be executed independently or as a complete workflow.
It performs genome alignments, CNE detection, evolutionary rate estimation.

---

**Global Options:**

| Option          | Description                           |
| --------------- | ------------------------------------- |
| `-h, --help`    | Show this help message and exit.      |
| `-v, --version` | Print the current version of CNEwrap. |

---

**Subcommands:**

| Subcommand | Description                                                                                                              |
| ---------- | ------------------------------------------------------------------------------------------------------------------------ |
| **align**  | Perform whole-genome alignment using **LASTZ**.                                                                          |
| **merge**  | Merge pairwise MAF alignments into **multi-way alignments** that include all species.                                    |
| **scne**   | Scan the entire aligned genome to identify **conserved noncoding elements (CNEs)** using **GERP** and **Phast** methods. |
| **trace**  | Extract and manipulate CNE alignments and sequence information for downstream analysis.                                  |
| **evolve** | Detect **accelerated CNEs** for specific species or clades.                                                              |
| **allrun** | Execute all steps sequentially — equivalent to running the full pipeline.                                                |

---


### 1. whole genome alignment with lastz

The `align` module performs **whole-genome alignment** among all species listed in the input tree.  
It uses **LASTZ** as the alignment engine to generate pairwise alignment files between the reference genome and all other genomes.  
Please note that the prefix of each FASTA file in GENOMEDIR must correspond exactly to the species name in the input tree. 
For local runs, -n controls the number of genome segments used for parallel processing.
Use --dryrun to verify the generated commands before launching large-scale alignment jobs.
For SGE runs, make sure --invoke_sge and --maxjobs are configured appropriately to match your cluster’s capacity.

<pre>
python $CNEHOME/cnewrap.py align -h 
usage: CNEwrap align [-h] [-r REFGENOME] [-t TREEFILE] [-g GENOMEDIR] [-n SPLITN] [-d {far,medium,near,lastz.dist}] [--invoke_sge] [--sge_para SGE_PARA]
                     [--maxjobs MAXJOBS] [--interval CHECK_INTERVAL] [--dryrun]

options:
  -h, --help            show this help message and exit
  -r REFGENOME          reference genome for aligning to (default: None)
  -t TREEFILE           treefile containing all species involved in analysis (default: None)
  -g GENOMEDIR          Directory containing all genome sequence files (default: None)
  -n SPLITN             parts of reference gnome to split into (default: 20)
  -d {far,medium,near,lastz.dist}
                        distance between species candidate: near medium far (default: medium)
  --invoke_sge          run lastz in batch with SGE cluster (default: False)
  --sge_para SGE_PARA   resource parameters that pass to SGE cluster qsub, valid only when SGE invoked (default: -cwd -l vf=20G,p=1 -pe smp 1)
  --maxjobs MAXJOBS     maxmaium number of concurrent jobs, valid only when SGE invoked (default: 100)
  --interval CHECK_INTERVAL
                        time interval to check status of concurrent jobs,valid only when SGE invoked (default: 10)
  --dryrun              runing without implementing it in real conditions (default: False)

</pre>

---

### 2. merge pairsise maf to multi-way alignments containing all species

The **`merge`** module is responsible for combining pairwise alignment results (generated by the `align` step) into multi-species multiple alignment files (MAFs).  
It automatically manages chromosome renaming, file linking, and incremental merging, ensuring that the resulting alignment includes all target species specified in the input phylogenetic tree. 
Ensure that the TREEFILE includes all species that were aligned in the previous step.
When using --renamechr, chromosome names will be permanently modified for compatibility across genomes. 
If merging was previously interrupted, using --increment allows resuming the process without restarting.

<pre>
python $CNEHOME/cnewrap.py merge -h 
usage: CNEwrap merge [-h] [-r REFGENOME] [-t TREEFILE] [-m MAFDIR] [-s SUFFIX] [-d TARGETDIR] [--renamechr] [--bypass_rename] [--increment] [--invoke_sge]
                     [--sge_para SGE_PARA] [--maxjobs MAXJOBS] [--check_interval CHECK_INTERVAL] [--dryrun]

options:
  -h, --help            show this help message and exit
  -r REFGENOME          species of reference genome for aligning to (default: None)
  -t TREEFILE           treefile containing all species involved in analysis (default: None)
  -m MAFDIR             merged maf by split parts (default: AllMAF)
  -s SUFFIX             suffix of final merged filtered maf files (default: net.filter.axt.maf)
  -d TARGETDIR          directory containning results of align lastz results, defaut: target (default: target)
  --renamechr           add species tag before chromosome and rename maffile together instead of softlink (default: False)
  --bypass_rename       skip step of rename or softlink maffiles (default: False)
  --increment           continue to rename or softlink maffiles instead of redo all maffiles (default: False)
  --invoke_sge          run lastz in batch with SGE cluster (default: False)
  --sge_para SGE_PARA   resource parameters that pass to SGE cluster qsub, valid only when SGE invoked (default: -cwd -l vf=1G,p=1)
  --maxjobs MAXJOBS     maxmaium number of concurrent jobs, valid only when SGE invoked (default: 100)
  --check_interval CHECK_INTERVAL
                        time interval to check status of concurrent jobs,valid only when SGE invoked (default: 30)
  --dryrun              runing without implementing it in real conditions (default: False)

</pre>

---

### 3. scan whole aligned gnome to identify CNEs with two method: GERP and Phast

The **`scne`** module is responsible for detecting **conserved noncoding elements (CNEs)** across multiple species based on the multi-species genome alignment results obtained in the previous steps.  
It integrates **GERP** and **PhastCons** (via the PHAST toolkit) to estimate sequence conservation and identify functionally constrained noncoding regions.  

If a precomputed mafinfo file is provided (--mafinfofile), the program will skip MAF parsing and directly identify CNEs.

<pre>
python $CNEHOME/cnewrap.py scne -h 
usage: CNEwrap scne [-h] [-r REFGENOME] [-t TREEFILE] [-m SPLITMAF] [-b MAFBLOCK] [-p PHASTDIR] [-g GERPDIR] [--phylopmod PHYLOPMOD] [--mafinfofile MAFINFOFILE]
                    [--method {GERP,Phast,both}] [--invoke_sge] [--sge_para SGE_PARA] [--maxjobs MAXJOBS] [--ncpu NCPU] [--check_interval CHECK_INTERVAL]

options:
  -h, --help            show this help message and exit
  -r REFGENOME          species of reference genome for aligning to (default: None)
  -t TREEFILE           treefile containing all species involved in analysis (default: None)
  -m SPLITMAF           seperated mafs one block per file (default: splitmaf)
  -b MAFBLOCK           directory of results of mafs one sequence alignment per file (default: MAFBlock)
  -p PHASTDIR           output of PhastCon results (default: phast_dir)
  -g GERPDIR            output of GERP results (default: gerp_dir)
  --phylopmod PHYLOPMOD
                        mod file produced by phyloFit (default: )
  --mafinfofile MAFINFOFILE
                        mafinfo file used for CNE identification, will override the -b and -m if setted (default: )
  --method {GERP,Phast,both}
                        method to identify CNEs (default: both)
  --invoke_sge          run lastz in batch with SGE cluster (default: False)
  --sge_para SGE_PARA   resource parameters that pass to SGE cluster qsub, valid only when SGE invoked (default: -cwd -l vf=1G,p=1)
  --maxjobs MAXJOBS     maxmaium number of concurrent jobs, valid only when SGE invoked (default: 100)
  --ncpu NCPU           maxmaium number of concurrent jobs, valid when SGE is not invoked (default: 100)
  --check_interval CHECK_INTERVAL
                        time interval to check status of concurrent jobs,valid only when SGE invoked (default: 5)

</pre>

---

### 4. extract and manipulate CNE alignments and sequence information

The **`trace`** module is used to extract, manipulate, and annotate **CNE alignments and sequence information** after CNE identification.  
It links conservation scores (from **GERP** and **PhastCons**) with reference genome annotations, allowing users to trace the genomic context of each conserved region and prepare sequence files.
- Extracts **aligned sequences** corresponding to identified CNEs from multi-species MAF files.  
- Integrates **GERP** and/or **PhastCons** conservation results for accurate CNE boundary definition.  
- Generates **BED** files for genomic coordinate visualization and **FASTA** files for sequence-based analyses.  
- Associates CNEs with annotated genomic features (e.g., CDS, intron, intergenic regions) using reference **GFF** files.  

<pre>
usage: CNEwrap trace [-h] [-r REFGENOME] [-a GFFFILE] [-e BEDKEY] [-g GERP_DIR] [-p PHAST_DIR] [-s SPLITMAF] [-f FASDIR]

options:
  -h, --help    show this help message and exit
  -r REFGENOME  species reference genome for aligning to (default: None)
  -a GFFFILE    annotation gff file of the reference genome (default: None)
  -e BEDKEY     keyword of bedfile for {keyword}{minlen}.bed format for all CNEs (default: cne)
  -g GERP_DIR   directory of GERP results (default: gerp_dir)
  -p PHAST_DIR  directory of PHAST results (default: phast_dir)
  -s SPLITMAF   directory containing MAF files used to extract CNE alignments (default: splitmaf)
  -f FASDIR     directory of output files for CNE alignments (default: bed_fasta)

</pre>

---

### 5. identify accelerated CNEs for specific species/clades

The **`evolve`** module is designed to detect **accelerated evolution in conserved noncoding elements (CNEs)** for specific species or clades.  
It compares the evolutionary rate of CNEs in foreground species against background lineages using phylogenetic models, allowing users to pinpoint regions under potential positive selection or lineage-specific acceleration.  
The Supports multiple methods for evolutionary rate assessment phyloP, EvoAcc and ForwardGenomics.
Useful scripts to run PhyloAcc were also provide and can be found in $CNHOME/scripts.


<pre>
python $CNEHOME/cnewrap.py evolve -h 
usage: CNEwrap evolve [-h] [-f FGSPECIES] [-t TREEFILE] [-d FASDIR] [-p PHYLOPMOD] [-F FRACTION] [-G GAPS] [--evo_method {FG,SP,PP,SPPP,all}] [--distfile DISTFILE]

options:
  -h, --help            show this help message and exit
  -f FGSPECIES          forground species to detect rapid evolving CNEs, multispecies should be sperated by comma (default: None)
  -t TREEFILE           treefile containing all species involved in analysis (default: None)
  -d FASDIR             input file including intrested genes(one geneid/line) (default: None)
  -p PHYLOPMOD          mod file produced by phyloFit (default: None)
  -F FRACTION           fraction of not altnative sites for a locus in alignment (default: 1)
  -G GAPS               fraction of gaps for a locus in alignment (default: 1)
  --evo_method {FG,SP,PP,SPPP,all}
                        method to asscess acc CNEs (default: SPPP)
  --distfile DISTFILE   DNA replacematrix file (default: None)

</pre>

---

### run all steps in one command

CNEwrap offers a modular yet fully automated framework for:

* Genome alignment
* Multi-species CNE identification
* Detection of accelerated evolution across lineages
* Extraction and analysis of conserved regions

Whether running step-by-step or via the **`allrun`** shortcut, CNEwrap provides an efficient and reproducible workflow for comparative genomics and evolutionary analysis.

<pre>
python ../cnewrap.py allrun -h 
usage: CNEwrap allrun [-h] -r REFGENOME [-g GENOMEDIR] -t TREEFILE -a GFFFILE [-f FGSPECIES] [-n SPLITN] [-p PHYLOPMOD] [-F FRACTION] [-G GAPS] [-s SUFFIX] [-d TARGETDIR]
                      [--divtime {far,medium,near,lastz.dist}] [--invoke_sge] [--sge_para SGE_PARA] [--maxjobs MAXJOBS] [--interval CHECK_INTERVAL] [--dryrun]
                      [--evo_method {FG,SP,PP,SPPP,all}] [--distfile DISTFILE]

options:
  -h, --help            show this help message and exit
  -r REFGENOME          reference species name (default: None)
  -g GENOMEDIR          Directory containing all genome sequence files (default: None)
  -t TREEFILE           species tree file (newick) (default: None)
  -a GFFFILE            reference GFF annotation file,keep only CDS items (default: None)
  -f FGSPECIES          forground species to detect rapid evolving CNEs, multispecies should be sperated by comma (default: None)
  -n SPLITN             parts of reference gnome to split into (default: 20)
  -p PHYLOPMOD          mod file produced by phyloFit (default: None)
  -F FRACTION           fraction of not altnative sites for a locus in alignment (default: 1)
  -G GAPS               fraction of gaps for a locus in alignment (default: 1)
  -s SUFFIX             suffix of final merged filtered maf files (default: net.filter.axt.maf)
  -d TARGETDIR          directory containning results of align lastz results, defaut: target (default: target)
  --divtime {far,medium,near,lastz.dist}
                        distance between species candidate: near medium far (default: medium)
  --invoke_sge          run lastz in batch with SGE cluster (default: False)
  --sge_para SGE_PARA   resource parameters that pass to SGE cluster qsub, valid only when SGE invoked (default: -cwd -l vf=20G,p=1 -pe smp 1)
  --maxjobs MAXJOBS     maxmaium number of concurrent jobs, valid only when SGE invoked (default: 100)
  --interval CHECK_INTERVAL
                        time interval to check status of concurrent jobs,valid only when SGE invoked (default: 10)
  --dryrun              runing without implementing it in real conditions (default: False)
  --evo_method {FG,SP,PP,SPPP,all}
                        method to asscess acc CNEs (default: SPPP)
  --distfile DISTFILE   DNA replacematrix file (default: None)

</pre>

---

## 💡 Example

Run the full pipeline with default parameters:

```bash
cnewrap allrun -r refsp -i genomes_dir -t CNE.tre -g refsp.cds.gff
```

<pre>
cd CNEWRAP/example/
./run_example.sh
</pre>
