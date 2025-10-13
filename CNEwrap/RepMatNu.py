#!/usr/bin/env python3
import sys
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio.Align import substitution_matrices as SubsMat
from Bio.Align import MultipleSeqAlignment

import numpy as np

def upper_align(c_align):
	return MultipleSeqAlignment([rec.upper() for rec in c_align])

def expand_matrix(original_matrix,target_labels=['A','C','G','T']):
	expanded_matrix = np.zeros((4, 4))
	original_labels=list(original_matrix.alphabet)
	#print(original_labels)
	for i, row_label in enumerate(original_labels):
		for j, col_label in enumerate(original_labels):
			if row_label  in target_labels:
				row_idx = target_labels.index(row_label)
			if col_label in target_labels:
				col_idx = target_labels.index(col_label)
			expanded_matrix[row_idx, col_idx] = original_matrix[i, j]
	return expanded_matrix
def replace_matrix(maffile,setcounts=1000,minlen=500):
	align_counts=0
	if maffile.endswith('maf'):
		fmt="maf"
	elif maffile.endswith('fas') or maffile.endswith('fa') or maffile.endswith('fasta'):
		fmt="fasta"
	else:
		print("unkown format of alignment files %s"%maffile,file=sys.stderr)
		sys.exit(3)
		
	c_aligns = AlignIO.parse(maffile, fmt)
	observed_frequencies=np.zeros((4,4))
	for c_align in c_aligns:
		if c_align.alignment.length >= minlen:
			c_align=upper_align(c_align)
			c_align_freq=c_align.substitutions
			if c_align_freq.shape != (4,4):
				c_align_freq=expand_matrix(c_align_freq)
			observed_frequencies  = observed_frequencies + c_align_freq 
			align_counts+=1
		if align_counts > setcounts:
			break
	print("%d alignments meet the criteria will be used to count replace matrix"%align_counts,file=sys.stderr)
	
	observed_frequencies /= np.sum(observed_frequencies)
	residue_frequencies = np.sum(observed_frequencies, 0)
	sum(residue_frequencies) == 1.0
	expected_frequencies = np.dot(residue_frequencies[:, None], residue_frequencies[None, :])
	my_lom = np.log2(observed_frequencies / expected_frequencies)
	return my_lom

if __name__ == "__main__":
	import sys
	maffile = "splitmaf/Tele161.maf"
	maffile=sys.argv[1]
	my_lom=replace_matrix(maffile)
	np.set_printoptions(precision=3) 
	print(my_lom)
	with open('dna.dist','w') as DistMX:
		print(my_lom.format("%6.3f"),file=DistMX)
