#!/usr/bin/env python3
import sys,os,argparse
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from glob import glob
from CNEwrap.utils import runpool
RESET = '\033[0m'
GREEN = '\033[1;32m'

def merge_alignments(alignment1, alignment2):
	# Check that both alignments have the same number of sequences
	if len(alignment1) != len(alignment2):
		raise ValueError("The number of sequences in both alignments must be the same.")

	# Concatenate each sequence, ensuring we maintain the annotations for 'start' and 'size'
	merged_records = []
	recNO=1
	#to count record in alignment, the first record should be the reference
	#and will be standard to concatenate alignment
	for record1, record2 in zip(alignment1, alignment2):
		# Merge the sequences by concatenating the sequences from both alignments
		merged_seq = record1.seq + record2.seq
		
		# Create a new SeqRecord for the merged sequence
		merged_record = record1[:]  # Copy the first record
		merged_record.seq = merged_seq
		
		# Update the annotations for 'start' and 'size'
		start1 = record1.annotations['start']
		srcSize1 = record1.annotations['srcSize']
		size1 = record1.annotations['size']
		strand1 = record1.annotations['strand']
		start2 = record2.annotations['start']
		srcSize2 = record2.annotations['srcSize']
		size2 = record2.annotations['size']
		strand2 = record2.annotations['strand']
		if recNO==1:
			assert strand1 == strand2, "strand of %s is not equal to concatenate"%record1.id
			assert srcSize1 == srcSize2, "srcSize of %s and %s is not equal to concatenate for %s and %s"%(srcSize1,srcSize2,record1,record2)
		elif strand1 != strand2 or srcSize1 != srcSize2:
			print("WARNING: strand or srcSize of are not equal to concatenate for %s and %s"%(record1,record2))

		merged_record.annotations['start'] = min(start1, start2)
		merged_record.annotations['size'] = size1+size2
		merged_record.annotations['strand'] = strand1
		merged_record.annotations['srcSize'] = srcSize1
		# Append the merged record to the list
		merged_records.append(merged_record)
		recNO+=1
	
	# Create a new MultipleSeqAlignment object from the merged records
	merged_alignment = MultipleSeqAlignment(merged_records)
	
	return merged_alignment

def reshape_maf(maffile,refsp,fmt="maf"):
	lastend=-1
	for alignment in AlignIO.parse(maffile,fmt):
		for record in alignment:
			chrom=record.id
			if chrom.startswith(refsp):
				chrlen=record.annotations['srcSize']
				start=record.annotations['start']
				size=record.annotations['size']
				strand=record.annotations['strand']
				assert strand == 1, "make sure refer sequence is plus strand"
				end=start+size
				#print(start,end,lastend)
				if start==lastend and chrom == lastchrom:
					alignment=merge_alignments(lastalign,alignment)
					#alignment=lastalign+alignment
					#print(alignment)
				else:
					if not lastend == -1:
						yield lastalign
						#print(lastalign)
				lastend=end
				lastchrom=chrom
		lastalign=alignment
	yield lastalign

def get_aln_pos(seqs,start,pos):
	#give sequence postion to calculate alignment position
	assert pos>=start, "gived position not in right interval, recheck your input"
	seqcursor=start
	alicursor=0
	while pos > seqcursor:
		alicursor+=1
		if seqs[alicursor]!="-":
				seqcursor+=1
	return alicursor #pos=index+1

def get_seq_pos_plus(seqs,start,pos):
	#give alignment postion to calculate sequence position
	assert pos>=0, "gived position not in right interval, recheck your input"
	seqcursor=0
	alicursor=0
	while alicursor < pos:
		alicursor+=1
		if seqs[alicursor]!="-":
			seqcursor+=1
	return start+seqcursor

def get_seq_pos_minus(seqs,end,pos):
	#give alignment postion to calculate sequence position
	assert pos>=0, "gived position not in right interval, recheck your input"
	seqcursor=0
	alicursor=0
	while alicursor < pos:
		alicursor+=1
		if seqs[alicursor]!="-":
			seqcursor+=1
	return end-seqcursor

def ExtractMafSeq(alignment,aimsp,refchr,refstart,refend):
	#get alignment start
	records={record.id.split("_")[0]:record for record in alignment}
	assert aimsp in records,"%s is not in alignment %s\n%s"%(aimsp,records,alignment)
	alignstart=get_aln_pos(records[aimsp].seq,records[aimsp].annotations['start'],refstart)
	alignend=get_aln_pos(records[aimsp].seq,records[aimsp].annotations['start'],refend+1)

	#get seqs start..end and sequences

	for record in alignment:
		chrom=record.id
		chrlen=record.annotations['srcSize']
		size=record.annotations['size']
		strand=record.annotations['strand']
		if strand == 1:
			start=record.annotations['start']
			end=start+size
			seqstart=get_seq_pos_plus(record.seq,start,alignstart)
			seqend=get_seq_pos_plus(record.seq,start,alignend)
		else:
			#ucscRevStart = chromSize - oneEnd
			#ucscRevEnd   = chromSize - (oneStart - 1)
			start=chrlen-record.annotations['start']-size-2
			end=chrlen-(record.annotations['start']-1)
			seqstart=get_seq_pos_minus(record.seq,end,alignstart)
			seqend=get_seq_pos_minus(record.seq,end,alignend)
		yield((record.id,seqstart,seqend,strand,alignstart,alignend,record.seq[alignstart:alignend]))

def process_bar(num, total):
	rate = float(num)/total
	ratenum = int(100*rate)
	r = '{}\r[{}{}]{}%{}'.format(GREEN,'#'*ratenum,' '*(100-ratenum), ratenum,RESET)
	sys.stdout.write(r)
	sys.stdout.flush()

def write_align(bedfile,maffile,aimsp,outdir,fmt="maf"):
	bedlist=[]
	with open(bedfile) as BED:
		for rec in BED:
			reclist=rec.strip().split()
			bedlist.append(reclist)
	if fmt == "maf":
		meltmaf=reshape_maf(maffile,aimsp,fmt)
		for alignment in meltmaf:
			for record in alignment:
				if record.id.startswith(aimsp):
					refchr=record.id
					#print(refchr,record)
					refstart=record.annotations['start']
					refend=record.annotations['start']+record.annotations['size']
				for bed in bedlist:
					if bed[0] == refchr:
						if int(bed[1]) >= refstart and int(bed[2]) <= refend-1:
							Fsign="G"+bed[3] if bed[3][0].isdigit() else bed[3]
							OUTFAS=open(os.sep.join([outdir,Fsign+".fas"]),'w')
							OUTGFF=open(os.sep.join([outdir,Fsign+".gff"]),'w')
							subalign=ExtractMafSeq(alignment,aimsp,refchr,int(bed[1]),int(bed[2])-1)
							for rec in subalign:
								#print(">%s\n%s"%(rec[0].split(".")[0],rec[6]),file=OUTFAS)
								print(">%s\n%s"%(rec[0].split("_")[0],rec[6]),file=OUTFAS)
								print("%s\t%s\t%s\t%s\t%s-%s"%rec[:-1],file=OUTGFF)
							OUTFAS.close
							OUTGFF.close
							bedlist.remove(bed)
	elif fmt == "fasta":
		for alignment in AlignIO.parse(maffile,fmt):
			for bed in bedlist:
				Fsign="G"+bed[3] if bed[3][0].isdigit() else bed[3]
				OUTFAS=open(os.sep.join([outdir.replace(".fasta","").replace(".fas",""),Fsign+".fas"]),'w')
				for record in alignment:
					subseq=record.seq[int(bed[1]):int(bed[2])]
					print(">%s\n%s"%(record.id.strip().split()[0],subseq),file=OUTFAS)
	else:
		sys.exit('format not recongnized')

def makeoutdir(maffile,fasdir):
	basemaf=os.path.basename(maffile)
	outdir=os.sep.join([fasdir,basemaf.replace(".melt.maf","").replace('.fasta','').replace('.fas','')])
	os.makedirs(outdir,exist_ok=True)
	return(outdir)

def batch_write_align(bedfile,mafdir,aimsp,fasdir="bed_fasta",fmt="maf"):
	maflist=glob(os.path.join(mafdir,"*.maf"))
	argslist=[(bedfile,maffile,aimsp,makeoutdir(maffile,fasdir),fmt) for maffile in maflist]
	runpool(write_align,argslist,len(argslist))
						
if __name__=="__main__":
#information from bedfile
	parser = argparse.ArgumentParser(description='Parameter configuration for genomic processing')
	parser.add_argument('--bedfile', type=str, default='cne30.bed',
		help='Default BED file path')
	parser.add_argument('--aimsp', type=str, default='simHuman',
		help='Species identifier for alignment')
	parser.add_argument('--fasdir', type=str, default='bed_fasta',
		help='Directory for FASTA files')
	parser.add_argument('--maffile', type=str, default="",
		help='MAF alignment file path')
	parser.add_argument('--mafdir', type=str, default='',
		help='Directory containing MAF files')
	args = parser.parse_args()

	#bedfile="cne30.bed"
	#aimsp="simHuman"
	#fasdir="bed_fasta"
	#maffile="splitmaf/simHuman1chr7.maf"
	#mafdir="splitmaf"
	#bedfile="/data/nfs2/lrh/Tbai_inputs/SimRat/Simulation_ratite/simu_500_200_diffr_2-6.bed3"
	#maffile="/data/nfs2/lrh/Tbai_inputs/SimRat/Simulation_ratite/simu_500_200_diffr_2-6.fasta"
	#aimsp="taeGut"
	bedfile=args.bedfile
	fasdir=args.fasdir
	aimsp=args.aimsp
	if args.maffile:
		maffile=args.maffile
		outdir=makeoutdir(maffile,fasdir)
		write_align(bedfile,maffile,aimsp,outdir,fmt="fasta")
	if args.mafdir:
		mafdir=args.mafdir
		batch_write_align(bedfile,mafdir,aimsp,fasdir,fmt="maf")
