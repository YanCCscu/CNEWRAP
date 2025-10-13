***
# CNEwrap A integrated pipeline for estimating Convergent accelerated evolutionary conserved non-coding elements
***
**CNEwrap** is a Python-based pipeline designed to identify conserved noncoding elements (CNEs) across multiple genomes, perform multiple sequence alignments, and conduct evolutionary analyses.
## Install
<pre>
git clone https://github.com/YanCCscu/CNEWRAP.git
conda env create -f cnewrap_env.yml
conda activate cnewrap_env
</pre>
Alternatively, you can install with conda command:
<pre>
conda install yccscucib::cnewrap
</pre>
## Useage
<pre>
python cnewrap.py -h
options:
  -h, --help            show this help message and exit
  -v, --version         print the version

subprocess:
  {align,merge,scne,trace,evolve,allrun}
    align               whole genome alignment with lastz
    merge               merge pairsise maf to multi-way alignments containing all species
    scne                scan whole aligned gnome to identify CNEs with two method: GERP and Phast
    trace               extract and manipulate CNE alignments and sequence information
    evolve              identify accelerated CNEs for specific species/clades
    allrun              run all steps
</pre>

### step 1: whole genome alignment with lastz
<pre>
python cnewrap.py align -h
usage: CNEwrap align [-h] [-r REFGENOME] [-t TREEFILE] [-g GENOMEDIR] [-n SPLITN] [-d {far,medium,near,lastz.dist}] [--invoke_sge] [--sge_para SGE_PARA] [--maxjobs MAXJOBS]
                     [--interval CHECK_INTERVAL] [--dryrun]

options:
  -h, --help            show this help message and exit
  -r REFGENOME          reference genome for aligning to (default: None)
  -t TREEFILE           treefile containing all species involved in analysis (default: None)
  -g GENOMEDIR          Directory containing all genome sequence files (default: None)
  -n SPLITN             parts of reference gnome to split into (default: 20)
  -d {far,medium,near,lastz.dist}
                        distance between species candidate: near medium far (default: medium)
  --invoke_sge          run lastz in batch with SGE cluster (default: False)
  --sge_para SGE_PARA   resource parameters that pass to SGE cluster qsub, valid only when SGE invoked (default: -cwd -l vf=40G,p=20 -pe smp 20)
  --maxjobs MAXJOBS     maxmaium number of concurrent jobs, valid only when SGE invoked (default: 100)
  --interval CHECK_INTERVAL
                        time interval to check status of concurrent jobs,valid only when SGE invoked (default: 10)
  --dryrun              runing without implementing it in real conditions (default: False)
</pre>


### step 2: merge pairsise maf to multi-way alignments containing all species

<pre>
python cnewrap.py merge -h
usage: CNEwrap merge [-h] [-r REFGENOME] [-t TREEFILE] [-m MAFDIR] [-s SUFFIX] [-d TARGETDIR] [--renamechr] [--bypass_rename] [--invoke_sge] [--sge_para SGE_PARA] [--maxjobs MAXJOBS]
                     [--check_interval CHECK_INTERVAL] [--dryrun]

options:
  -h, --help            show this help message and exit
  -r REFGENOME          species of reference genome for aligning to (default: None)
  -t TREEFILE           treefile containing all species involved in analysis (default: None)
  -m MAFDIR             merged maf by split parts (default: AllMAF)
  -s SUFFIX             suffix of final merged filtered maf files (default: net.filter.axt.maf)
  -d TARGETDIR          directory containning results of align lastz results, defaut: target (default: target)
  --renamechr           add species tag before chromosome and rename maffile together instead of softlink (default: False)
  --bypass_rename       skip step of rename or softlink maffiles (default: False)
  --invoke_sge          run lastz in batch with SGE cluster (default: False)
  --sge_para SGE_PARA   resource parameters that pass to SGE cluster qsub, valid only when SGE invoked (default: -cwd -l vf=1G,p=1)
  --maxjobs MAXJOBS     maxmaium number of concurrent jobs, valid only when SGE invoked (default: 100)
  --check_interval CHECK_INTERVAL
                        time interval to check status of concurrent jobs,valid only when SGE invoked (default: 30)
  --dryrun              runing without implementing it in real conditions (default: False)
</pre>

### step 3: scan whole aligned gnome to identify CNEs with two method: GERP and Phast

<pre>
python cnewrap.py scne -h
usage: CNEwrap scne [-h] [-r REFGENOME] [-t TREEFILE] [-m SPLITMAF] [-b MAFBLOCK] [-p PHASTDIR] [-g GERPDIR] [--method {GERP,Phast,both}] [--invoke_sge] [--sge_para SGE_PARA]
                    [--maxjobs MAXJOBS] [--check_interval CHECK_INTERVAL]

options:
  -h, --help            show this help message and exit
  -r REFGENOME          species of reference genome for aligning to (default: None)
  -t TREEFILE           treefile containing all species involved in analysis (default: None)
  -m SPLITMAF           seperated mafs one block per file (default: splitmaf)
  -b MAFBLOCK           directory of results of mafs one sequence alignment per file (default: MAFBlock)
  -p PHASTDIR           output of PhastCon results (default: phast_dir)
  -g GERPDIR            output of GERP results (default: gerp_dir)
  --method {GERP,Phast,both}
                        method to identify CNEs (default: both)
  --invoke_sge          run lastz in batch with SGE cluster (default: False)
  --sge_para SGE_PARA   resource parameters that pass to SGE cluster qsub, valid only when SGE invoked (default: -cwd -l vf=1G,p=1)
  --maxjobs MAXJOBS     maxmaium number of concurrent jobs, valid only when SGE invoked (default: 100)
  --check_interval CHECK_INTERVAL
                        time interval to check status of concurrent jobs,valid only when SGE invoked (default: 5)
</pre>
### step 4: extract and manipulate CNE alignments and sequence information
<pre>
python cnewrap.py trace -h
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

### step 5: identify accelerated CNEs for specific species/clades
<pre>
python cnewrap.py evolve -h
usage: CNEwrap evolve [-h] [-f FGSPECIES] [-t TREEFILE] [-d FASDIR] [-p PHYLOPMOD] [-F FRACTION] [-G GAPS] [--evo_method {FG,SP,PP,all}] [--distfile DISTFILE]

options:
  -h, --help            show this help message and exit
  -f FGSPECIES          forground species to detect rapid evolving CNEs, multispecies should be sperated by comma (default: None)
  -t TREEFILE           treefile containing all species involved in analysis (default: None)
  -d FASDIR             input file including intrested genes(one geneid/line) (default: None)
  -p PHYLOPMOD          mod file produced by phyloFit (default: None)
  -F FRACTION           fraction of not altnative sites for a locus in alignment (default: 1)
  -G GAPS               fraction of gaps for a locus in alignment (default: 1)
  --evo_method {FG,SP,PP,all}
                        method to asscess acc CNEs (default: all)
  --distfile DISTFILE   path of DNA distance file (default: None)

</pre>

### example of running CNEWRAP
<pre>
cd CNEWRAP/example/
sh run_example.sh
</pre>
