#!/usr/bin/env python3
import sys
import itertools
from Bio import AlignIO
from math import log
from glob import glob
from collections import Counter
from CNEwrap.utils import runpool
from ete3 import Tree
SIG_THRESHOLD=4

def calc_group_distance_normalized(treefile, species):
	t = Tree(treefile, format=1)
	all_species = [leaf.name for leaf in t.get_leaves()]
	foreground_species=species.split(",")
	fg_species = list(set(foreground_species) & set(all_species))
	bg_species = list(set(all_species) - set(fg_species))

	if len(fg_species) == 0 or len(bg_species) == 0:
		raise ValueError("Error, no forground species offered,")

	total_dist = 0
	count = 0
	for fg in fg_species:
		for bg in bg_species:
			dist = t.get_distance(fg, bg)
			total_dist += dist
			count += 1
	avg_fg_bg = total_dist / count if count > 0 else 0
	max_dist = 0
	for sp1, sp2 in itertools.combinations(all_species, 2):
		dist = t.get_distance(sp1, sp2)
		if dist > max_dist:
			max_dist = dist

	norm_score = avg_fg_bg / max_dist if max_dist > 0 else 0
	return avg_fg_bg, norm_score, max_dist

def information_content(seqcol): #Thomas et al., 1990
	uniAA=set(seqcol)
	lenAA=len(seqcol)
	H=0
	for AA in uniAA:
		propAA=seqcol.count(AA)/lenAA
		H+=propAA*log(propAA,2)
	return(2-(-H))

#tell whethe the position is species spcific
def sp_acceleration_score(seqcol, aimAA, outAA):
	total_length = len(seqcol)
	aim_counter = Counter(aimAA)
	out_counter = Counter(outAA)
	aim_set = set(aim_counter.keys())
	out_set = set(out_counter.keys())
	aim_consistency = 0
	#Jaccard distance 
	intersection = aim_set & out_set
	union = aim_set | out_set
	if union:
		# Jaccard distance = 1 - (intersection/ union)
		diff_score = 1 - len(intersection) / len(union)
	aim_consistency=information_content(aimAA)
	out_consistency=information_content(outAA)
	score = diff_score*(aim_consistency * (len(aimAA)/total_length))*(out_consistency * (len(outAA)/total_length)) 
	return round(score, 4)

def weighted_AA_distance(aimAA, outAA, DNAorder, distdat):
	total_weighted_dist = 0
	total_weight = 0
	aim_counter=Counter(aimAA)
	out_counter=Counter(outAA)
	for aim_aa, aim_freq in aim_counter.items():
		for out_aa, out_freq in out_counter.items():
			distance = DNAdist(aim_aa, out_aa, DNAorder, distdat)
			weight = (aim_freq / sum(aim_counter.values())) * (out_freq / sum(out_counter.values()))
			total_weighted_dist += distance * weight
			total_weight += weight
	if total_weight == 0:
		return -6
	return total_weighted_dist / total_weight

def sp_specific(seqcol,aimAA,outAA,fraction=1,gapmax=1):
	maxaimAA=max(aimAA,key=aimAA.count)
	maxoutAA=max(outAA,key=outAA.count)
	if gapmax>=1:
		if seqcol.count('-')>gapmax:
			return False
	elif 1>gapmax>0:
		if seqcol.count('-')/len(seqcol)>gapmax:
			return False
	else:
		raise Exception("gapmax should be a number >0")
			
	if len(set(aimAA)) >= 3:
		return False
	if aimAA.count(maxaimAA)>=len(aimAA)*float(fraction) and not maxaimAA in outAA:
		return True
	#leftAA=seqcol.replace(aimAA[0],'')
	for aa in set(aimAA):
		if aa != "-" and aa in outAA:
			return False
	return True
def flankAA(alignment,splist,sp_dict,i,f=10):
	fAAset=[]
	alignment_length=alignment.get_alignment_length()
	fstart = i-f if i-f >= 0 else 0
	fend = i+f+1 if i+f+1 < alignment_length else alignment_length
	for ii in range(fstart,fend):
		seqcol=alignment[:,ii]
		fAA=[seqcol[sp_dict[ai]] for ai in splist]
		fAAset+=fAA
	return(fAAset)
def parse_dist(distfile):																					
	distlist=[]																					  
	with open(distfile) as DF:																			   
		header=DF.readline()																			 
		hlist=header.strip().split()																		 
		for line in DF:																				  
			llist=line.strip().split()																	   
			distlist.append([float(l) for l in llist[1:]])															   
	return (hlist,distlist)						 

def DNAdist(A1,A2,DNAorder,distdat):
	DNAindex={a:i for i,a in enumerate(DNAorder) }
	if A1 in DNAindex and A2 in DNAindex:
		return(distdat[DNAindex[A1]][DNAindex[A2]])
	else:
		return(-6)

def SpSpSite(alifile,species,treefile,DNAorder,distdat,bgfile=None,SP=sys.stdout,SPSUM=sys.stdout,fraction=1,gaps=1):
	spsp_scores=[]
	alignment = AlignIO.read(alifile, "fasta")
	alignment_length=alignment.get_alignment_length()
	sp_order=[record.id.split('|')[0] for record in alignment]
	sp_dict={sp:i for i,sp in enumerate(sp_order)}
	splist=species.strip().split(",")
	if (set(splist) & set(sp_order)) != set(splist):
		print("PASS: %s not all specific species found in alignment"%alignment,end="\t")
		splist= list(set(splist) & set(sp_order))
		print("will use %s as forground species"%splist)
		#sys.exit(2)
	outlist=list(set(sp_order)-set(splist))
	if bgfile:
		with open(bgfile) as BG:
			bgsps=[sp.strip() for sp in BG]
			outlist=list((set(sp_order) & set(bgsps))-set(splist))
	raw_dist, norm_dist, max_dist = calc_group_distance_normalized(treefile, species)
	#print(f" : {norm_dist:.4f}")
	for i in range(alignment_length):
		seqcol=alignment[:,i]
		aimAA=[seqcol[sp_dict[ai]].upper() for ai in splist if ai in sp_order]
		if '-' in aimAA:
			fAAset=flankAA(alignment,splist,sp_dict,i,f=10)
			#print(fAAset,file=sys.stderr)
			if fAAset.count('-')/len(fAAset)>0.5 or fAAset.count('-')/(alignment_length*(10+1))>len(splist)/alignment_length:
				continue
		outAA=[seqcol[sp_dict[oi]].upper() for oi in outlist]
		acc_score=sp_acceleration_score(seqcol,aimAA,outAA)
		if acc_score > 0:
			siteIC=information_content(seqcol)
			aadist=weighted_AA_distance(aimAA, outAA, DNAorder, distdat)
			sigvalue=(siteIC+aadist/(-3)+acc_score)/norm_dist
			spsp_scores.append(sigvalue)
			usedsps=[ai for ai in splist if ai in sp_order]
			print("%s\t%s\t%.3f\t%.3f\t%.3f\t%.3f\t%s\t%s\t%s"%(alifile,i+1,siteIC,aadist,acc_score,sigvalue,\
			"".join(sorted(list(set(outAA)),reverse=True)),set(aimAA),",".join(usedsps)),file=SP)
	#get alignment information
	varsites=len(spsp_scores)
	if varsites:
		spsp_score_sum=sum(spsp_scores)
		spsp_score_max=max(spsp_scores)
		spsp_score_ave=spsp_score_sum/varsites
		signal_acc="ACC" if spsp_score_max > SIG_THRESHOLD or spsp_score_ave > 2.5 else "NUT"
		print("{alifile}\t{seqlen}\t{varsites}\t{spsp_score}\t{max_score}\t{ave_score}\t{acc}".format(alifile=alifile,
				seqlen=alignment_length,
				varsites=varsites,
				spsp_score=spsp_score_sum,	
			max_score=spsp_score_max,
			ave_score=round(spsp_score_ave,4),
			acc=signal_acc),file=SPSUM)

def batch_SpSpSite(alidir,species,treefile,DNAorder,distdat,bgfile=None,SP=sys.stdout,SPSUM=sys.stdout,fraction=1,gaps=1):
	alifilelist=glob("%s/*.fas"%alidir)
	argslist=[(alifile,species,DNAorder,distdat,bgfile,SP,SPSUM,fraction,gaps) for alifile in alifilelist if "anc.fas" not in alifile]
	ncpu=len(argslist) if len(argslist)<20 else 20
	runpool(SpSpSite,argslist,ncpu)

def batch_SpSpSite2(alidir,species,treefile,distfile,bgfile=None,SP=sys.stdout,SPSUM=sys.stdout,fraction=1,gaps=1):
	alifilelist=glob("%s/*.fas"%alidir)
	(DNAorder,distdat)=(parse_dist(distfile))
	for alifile in alifilelist:
		SpSpSite(alifile,species,treefile,DNAorder,distdat,bgfile=None,SP=SP,SPSUM=SPSUM,fraction=1,gaps=1)

if __name__ == "__main__":
	alidir="bed_fasta/simu_AccelerationPattern_CaseA"
	species="aptHaa,aptOwe,aptRow,casCas,droNov"
	distfile="dna.dist"
	treefile="./ratite.tre"
	OUTFILE=open("spout.txt","w")
	SUMFILE=open("spsum.txt","w")
	batch_SpSpSite2(alidir,species,treefile,distfile,bgfile=None,SP=OUTFILE,SPSUM=SUMFILE,fraction=1,gaps=1)

