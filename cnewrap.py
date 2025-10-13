#!/usr/bin/env python3
# A pipeline for dealing with CNE discovery and evolutionary charactor
import argparse,os,sys
from CNEwrap.align import genome_align
from CNEwrap.merge import merge_maf
from CNEwrap.scan import cne_scan
from CNEwrap.trace import trace_seqs 
from CNEwrap.evolve import acc_cne
VERSION="1.2"
cmdir=os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])),"bin")
os.environ['DRMAA_LIBRARY_PATH']=os.path.join(cmdir,"lib","libdrmaa.so.1.0")
os.system('export DRMAA_LIBRARY_PATH')
if __name__ == '__main__':
	#program_dir=os.path.dirname(os.path.abspath(sys.argv[0]))
#############parse arguments#############
	parser = argparse.ArgumentParser(prog="CNEwrap",epilog="example: cnewrap allrun -r abc -i genomes",
	description='CNEWRAP: an integrated pipeline for estimating accelerated evolutionary conserved non-coding elements')
	subparsers = parser.add_subparsers(title="subprocess",dest='command')
	parser.add_argument("-v", "--version", action="store_true",dest="version", help="print the version")
	#parser.add_argument("-r", "--refgenome", action="store",dest="refgenome", help="reference genome for aligning to")
	#parser.add_argument('-i', "--indir", action='store', dest="indir", help='Directory containing all genome sequence files for input')

#add main subprocesses
	align = subparsers.add_parser('align', help='whole genome alignment with lastz',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	merge = subparsers.add_parser('merge', help='merge pairsise maf to multi-way alignments containing all species',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	scne = subparsers.add_parser('scne', help='scan whole aligned gnome to identify CNEs with two method: GERP and Phast',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	trace = subparsers.add_parser('trace', help='extract and manipulate CNE alignments and sequence information',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	evolve = subparsers.add_parser('evolve', help='identify accelerated CNEs for specific species/clades',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	allrun = subparsers.add_parser('allrun', help='run all steps',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
#command 'align'
	"""
	genome_align(genomedir,refgenome,splitn,divtime,cmdir)
	"""
	align.add_argument('-r', action='store', dest="refgenome", help='reference genome for aligning to')
	align.add_argument('-t', action='store', dest="treefile", help='treefile containing all species involved in analysis')
	align.add_argument('-g', action='store', dest="genomedir", help='Directory containing all genome sequence files')
	align.add_argument('-n', action='store', dest="splitn", default=20,help='parts of reference gnome to split into')
	align_choices=["far","medium","near","lastz.dist"]
	align.add_argument('-d', action='store', dest="divtime", choices=align_choices,default='medium', help='distance between species candidate: near medium far')
	align.add_argument('--invoke_sge', action='store_true', help='run lastz in batch with SGE cluster')
	align.add_argument('--sge_para', action='store', dest="sge_para",default="-cwd -l vf=20G,p=1 -pe smp 1", \
		help='resource parameters that pass to SGE cluster qsub, valid only when SGE invoked')
	align.add_argument('--maxjobs', action='store', dest="maxjobs", default=100, help='maxmaium number of concurrent jobs, valid only when SGE invoked')
	align.add_argument('--interval', action='store', dest="check_interval", default=10, help='time interval to check status of concurrent jobs,valid only when SGE invoked')
	align.add_argument('--dryrun', action='store_true', help='runing without implementing it in real conditions')
	align.set_defaults(func=genome_align)

#command 'merge'
	"""
	merge_maf(treefile,refsp,suffix,AllMAF,cmdir):
	"""
	merge.add_argument('-r', action='store', dest="refgenome", help='species of reference genome for aligning to')
	merge.add_argument('-t', action='store', dest="treefile", help='treefile containing all species involved in analysis')
	merge.add_argument('-m', action='store', dest="mafdir", default='AllMAF',help='merged maf by split parts')
	merge.add_argument('-s', action='store', dest="suffix", default='net.filter.axt.maf',help='suffix of final merged filtered maf files')
	merge.add_argument('-d', action='store', dest="targetdir", default='target',help='directory containning results of align lastz results, defaut: target')
	merge.add_argument('--renamechr', action='store_true', help='add species tag before chromosome and rename maffile together instead of softlink')
	merge.add_argument('--bypass_rename', action='store_true', help='skip step of rename or softlink maffiles')
	merge.add_argument('--increment', action='store_true', help='continue to rename or softlink maffiles instead of redo all maffiles')
	merge.add_argument('--invoke_sge', action='store_true', help='run lastz in batch with SGE cluster')
	merge.add_argument('--sge_para', action='store', dest="sge_para", default="-cwd -l vf=1G,p=1", \
		help='resource parameters that pass to SGE cluster qsub, valid only when SGE invoked')
	merge.add_argument('--maxjobs', action='store', dest="maxjobs", default=100, help='maxmaium number of concurrent jobs, valid only when SGE invoked')
	merge.add_argument('--check_interval', action='store', dest="check_interval", default=30, help='time interval to check status of concurrent jobs,valid only when SGE invoked')
	merge.add_argument('--dryrun', action='store_true', help='runing without implementing it in real conditions')
	merge.set_defaults(func=merge_maf)

# command 'scne'
	"""
	cne_scan(mafblock,treefile,refgenome,phastdir,cmdir)
	"""
	scne.add_argument('-r', action='store', dest="refgenome", help='species of reference genome for aligning to')
	scne.add_argument('-t', action='store', dest="treefile", help='treefile containing all species involved in analysis')
	scne.add_argument('-m', action='store', dest="splitmaf", default='splitmaf',help='seperated mafs one block per file')
	scne.add_argument('-b', action='store', dest="mafblock",default="MAFBlock", help='directory of results of mafs one sequence alignment per file')
	scne.add_argument('-p', action='store', dest="phastdir",  default="phast_dir",help='output of PhastCon results')
	scne.add_argument('-g', action='store', dest="gerpdir",  default="gerp_dir",help='output of GERP results')
	scne.add_argument('--phylopmod', action='store', dest="phylopmod",default="", help='mod file produced by phyloFit')
	scne.add_argument('--mafinfofile', action='store', dest="mafinfofile",default="", help='mafinfo file used for CNE identification, will override the -b and -m if setted')
	scne.add_argument('--method', action='store', dest="cne_method",  default="both",choices=["GERP","Phast","both"],help="method to identify CNEs")
	scne.add_argument('--invoke_sge', action='store_true', help='run lastz in batch with SGE cluster')
	scne.add_argument('--sge_para', action='store', dest="sge_para", default="-cwd -l vf=1G,p=1", \
		help='resource parameters that pass to SGE cluster qsub, valid only when SGE invoked')
	scne.add_argument('--maxjobs', action='store', dest="maxjobs", default=100, help='maxmaium number of concurrent jobs, valid only when SGE invoked')
	scne.add_argument('--ncpu', action='store', dest="ncpu", default=100, type=int, help='maxmaium number of concurrent jobs, valid when SGE is not invoked')
	scne.add_argument('--check_interval', action='store', dest="check_interval", default=5, help='time interval to check status of concurrent jobs,valid only when SGE invoked')
	
	scne.set_defaults(func=cne_scan)

#command 'trace'
	"""
	trace_seqs(gerp_dir,phast_dir,splitmaf,gff_file,refgenome,bedfile="cne30.bed",fasdir="bed_fasta")
	"""
	trace.add_argument('-r', action='store', dest="refgenome", help='species reference genome for aligning to')
	trace.add_argument('-a', action='store', dest="gfffile", help='annotation gff file of the reference genome')
	trace.add_argument('-e', action='store', dest="bedkey", default="cne",  help='keyword of bedfile for {keyword}{minlen}.bed format for all CNEs')
	trace.add_argument('-g', action='store', dest="gerp_dir", default="gerp_dir",help='directory of GERP results')
	trace.add_argument('-p', action='store', dest="phast_dir",default="phast_dir", help='directory of PHAST results')
	trace.add_argument('-s', action='store', dest="splitmaf",default="splitmaf", help='directory containing MAF files used to extract CNE alignments')
	trace.add_argument('-f', action='store', dest="fasdir", default='bed_fasta', help='directory of output files for CNE alignments')
	trace.set_defaults(func=trace_seqs)

#command 'evolve'
	"""
	acc_cne(fasdir,mytree,refgenome,cmdir)
	"""
	evolve.add_argument('-f', action='store', dest="fgspecies", help='forground species to detect rapid evolving CNEs, multispecies should be sperated by comma')
	evolve.add_argument('-t', action='store', dest="treefile", help='treefile containing all species involved in analysis')
	evolve.add_argument('-d', action='store', dest="fasdir",  help='input file including intrested genes(one geneid/line)')
	evolve.add_argument('-p', action='store', dest="phylopmod", help='mod file produced by phyloFit')
	evolve.add_argument('-F', action='store', dest="fraction",default=1, help='fraction of not altnative sites for a locus in alignment')
	evolve.add_argument('-G', action='store', dest="gaps",default=1, help='fraction of gaps for a locus in alignment')
	evolve.add_argument('--evo_method', action='store', dest="evo_method",default="SPPP",choices=["FG","SP","PP","SPPP","all"], help='method to asscess acc CNEs')
	evolve.add_argument('--distfile', action='store', dest="distfile",default=None, help='DNA replacematrix file')
	evolve.set_defaults(func=acc_cne)

#define allbyone

	def allbyone(genomedir,refgenome,treefile,gfffile,fgspecies,splitn,divtime,cmdir,suffix,\
		mafdir="AllMAF",mafblock="MAFBlock",splitmaf="splitmaf",phastdir="phast_dir",gerpdir="gerp_dir",\
		fasdir="bed_fasta",bedkey="cne",\
		bgfile=None,fraction=1,gaps=1,evo_method="all",distfile=None,\
		invoke_sge=False,maxjobs=100,check_interval=10,dryrun=False):
		# Step 1. Alignment
		print("[1/5] Genome alignment ...")
		genome_align(genomedir,refgenome,treefile,splitn,divtime,cmdir,\
		invoke_sge=invoke_sge,sge_para="-l vf=4G,p=2 -pe smp 2",maxjobs=maxjobs,check_interval=check_interval,dryrun=dryrun)
		
		# Step 2. Merge maf
		print("[2/5] Merge MAF files ...")
		merge_maf(treefile,refgenome,suffix,cmdir,mafdir="AllMAF",targetdir="target",renamechr=False, \
		invoke_sge=invoke_sge,sge_para="-l vf=10G,p=1",maxjobs=maxjobs,check_interval=check_interval,\
		bypass_rename=False,increment=False, dryrun=dryrun)
		
		# Step 3. Scan for CNEs
		print("[3/5] Scan for CNEs ...")
		cne_scan(mafblock,treefile,refgenome,splitmaf,phastdir,gerpdir,cmdir,\
		mafinfofile="",phylopmod="",invoke_sge=invoke_sge,ncpu=40,maxjobs=maxjobs,check_interval=1,\
		sge_para="-l vf=1G,p=1",cne_method="both")
		
		# Step 4. Trace sequences
		print("[4/5] Trace CNE sequences ...")
		trace_seqs(gerpdir,phastdir,splitmaf,gfffile,refgenome,bedkey="cne",fasdir="bed_fasta")
		
		# Step 5. Evolutionary acceleration
		print("[5/5] Detect accelerated CNEs ...")
		acc_cne(fasdir,treefile,fgspecies,phylopmod,cmdir,bgfile=None,fraction=1,gaps=1,evo_method="all",distfile=None)
		print("=== [CNEwrap allrun finished successfully] ===")
#command 'TF bind motif'
	allrun.add_argument('-r', action='store', dest="refgenome", required=True, help='reference species name')
	allrun.add_argument('-g', action='store', dest="genomedir", help='Directory containing all genome sequence files')
	allrun.add_argument('-t', action='store', dest="treefile", required=True, help='species tree file (newick)')
	allrun.add_argument('-a', action='store', dest="gfffile", required=True, help='reference GFF annotation file,keep only CDS items')
	allrun.add_argument('-f', action='store', dest="fgspecies", help='forground species to detect rapid evolving CNEs, multispecies should be sperated by comma')
	allrun.add_argument('-n', action='store', dest="splitn", default=20,help='parts of reference gnome to split into')
	allrun.add_argument('-p', action='store', dest="phylopmod", help='mod file produced by phyloFit')
	allrun.add_argument('-F', action='store', dest="fraction",default=1, help='fraction of not altnative sites for a locus in alignment')
	allrun.add_argument('-G', action='store', dest="gaps",default=1, help='fraction of gaps for a locus in alignment')
	allrun.add_argument('-s', action='store', dest="suffix", default='net.filter.axt.maf',help='suffix of final merged filtered maf files')
	allrun.add_argument('-d', action='store', dest="targetdir", default='target',help='directory containning results of align lastz results, defaut: target')
	allrun.add_argument('--divtime', action='store', dest="divtime", choices=align_choices,default='medium', help='distance between species candidate: near medium far')
	allrun.add_argument('--invoke_sge', action='store_true', help='run lastz in batch with SGE cluster')
	allrun.add_argument('--sge_para', action='store', dest="sge_para",default="-cwd -l vf=20G,p=1 -pe smp 1", \
	help='resource parameters that pass to SGE cluster qsub, valid only when SGE invoked')
	allrun.add_argument('--maxjobs', action='store', dest="maxjobs", default=100, help='maxmaium number of concurrent jobs, valid only when SGE invoked')
	allrun.add_argument('--interval', action='store', dest="check_interval", default=10, help='time interval to check status of concurrent jobs,valid only when SGE invoked')
	allrun.add_argument('--dryrun', action='store_true', help='runing without implementing it in real conditions')
	allrun.add_argument('--evo_method', action='store', dest="evo_method",default="SPPP",choices=["FG","SP","PP","SPPP","all"], help='method to asscess acc CNEs')
	allrun.add_argument('--distfile', action='store', dest="distfile",default=None, help='DNA replacematrix file')
	allrun.set_defaults(func=allbyone)
	
	args = parser.parse_args() #deal with input 
	if not args.command:
		parser.print_usage()
		sys.exit(1)
	if args.version:
		print("CNEwrap v%s"%VERSION)
		sys.exit(1)
	if not hasattr(args, 'func'):
		args = parser.parse_args(['-h'])
	args_dict = {i: j for i, j in vars(args).items()}
	if args_dict['func'].__name__ != "trace_seqs":
		args_dict['cmdir']=cmdir
	args_dict.pop('func')
	args_dict.pop('version')
	args_dict.pop('command')
	if  not ('refgenome' in args_dict  or 'fgspecies' in args_dict):
		subparsers.choices[args.command].print_help()
	else:
		args.func(**args_dict)
