#!/usr/bin/env python3
import sys,os
import itertools
from Bio import AlignIO
from math import log
from glob import glob
from collections import Counter
from ete3 import Tree
import numpy as np
from scipy.stats import gamma
SIG_THRESHOLD=0.05
K_dist=1

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

def sp_acceleration_score(seqcol, FGBP, BGBP):
	total_length = len(seqcol)
	aim_counter = Counter(FGBP)
	out_counter = Counter(BGBP)
	aim_set = set(aim_counter.keys())
	out_set = set(out_counter.keys())
	aim_consistency = 0
	#Jaccard distance 
	intersection = aim_set & out_set
	union = aim_set | out_set
	if union:
		# Jaccard distance = 1 - (intersection/ union)
		diff_score = 1 - len(intersection) / len(union)
	fgIC=information_content(FGBP)
	bgIC=information_content(BGBP)
	siteIC=information_content(seqcol)
	score = diff_score*(fgIC)/(fgIC+bgIC)#(siteIC**2) #*(len(FGBP)/total_length)
	return (round(score, 4),diff_score,fgIC,bgIC,siteIC)

"""
def weighted_distance(FGBP, BGBP, DNAorder, distdat, gap_penalty=-6):
	#convert distda to dict
	score = {}
	for i, bi in enumerate(DNAorder):
		score[bi] = {}
		for j, bj in enumerate(DNAorder):
			#score[bi][bj] = distdat[i][j]
			score[bi][bj]=distdat[i][j]
	n_fg = len(FGBP)
	bg_counts = Counter(BGBP)
	N_bg = len(BGBP)
	S_bar = 0.0	  # ave score of background
	S_self = 0.0	 # self mutate score
	for b, cnt in bg_counts.items():
		#if b not in score:
		#	print("BGBP %s not exist in DNAorder" % b)
		p_b = float(cnt) / N_bg
		mean_b = sum(score.get(b,{}).get(x,gap_penalty) for x in FGBP) / n_fg
		S_bar += p_b * mean_b
		S_self += p_b * score.get(b, {}).get(b, gap_penalty)
	D = S_self - S_bar
	return {'S_bar': S_bar, 'S_self': S_self, 'D': D}
"""
def weighted_distance(
		FGBP, BGBP, DNAorder, distdat,
		gap_penalty=-6,		# gap to norm AA 
		gap_gap_score=3,	   # gap-gap
		n_match_score=3,	   # N-N 
		n_mismatch_score=-3	# N-other AA
	):
	# convert distdat to dict  
	score = {}
	for i, bi in enumerate(DNAorder):
		score[bi] = {}
		for j, bj in enumerate(DNAorder):
			score[bi][bj] = distdat[i][j]

	n_fg = len(FGBP)
	bg_counts = Counter(BGBP)
	N_bg = len(BGBP)

	def get_score(a, b):
		if a == '-' and b == '-':
			return gap_gap_score
		if a == '-' or b == '-':
			return gap_penalty
		if str(a).upper() == 'N' and str(b).upper() == 'N':
			return n_match_score
		if str(a).upper() == 'N' or str(b).upper() == 'N':
			return n_mismatch_score
		return score.get(a, {}).get(b, gap_penalty)

	S_bar = 0.0
	S_self = 0.0
	for b, cnt in bg_counts.items():
		p_b = float(cnt) / N_bg
		mean_b = sum(get_score(b, x) for x in FGBP) / n_fg
		S_bar += p_b * mean_b
		S_self += p_b * get_score(b, b)

	D = S_self - S_bar
	return {'S_bar': S_bar, 'S_self': S_self, 'D': D}


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

def SpSpSite(alifile,species,norm_dist,DNAorder,distdat,bgfile=None,fraction=1,gaps=1):
	spsp_scores=[]
	spsp_sites=[]
	spsp_paras=[]
	gid=os.path.basename(alifile).replace(".fas","").replace(".fasta","")
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
	#print(f" : {norm_dist:.4f}")
	for i in range(alignment_length):
		seqcol=alignment[:,i]
		FGBP=[seqcol[sp_dict[ai]].upper() for ai in splist if ai in sp_order]
		if '-' in FGBP:
			fAAset=flankAA(alignment,splist,sp_dict,i,f=10)
			#print(fAAset,file=sys.stderr)
			if fAAset.count('-')/len(fAAset)>0.5 or fAAset.count('-')/(alignment_length*(10+1))>len(splist)/alignment_length:
				continue
		BGBP=[seqcol[sp_dict[oi]].upper() for oi in outlist]
		(acc_score,diff_score,fgIC,bgIC,siteIC)=sp_acceleration_score(seqcol,FGBP,BGBP)
		if acc_score > 0:
			weightD=weighted_distance(FGBP, BGBP, DNAorder, distdat,gap_penalty=-6)
			dnadist=weightD["D"]
			##################
			sigvalue=(dnadist/(K_dist)*acc_score)/norm_dist
			##################
			#sigvalue=(K_dist/dnadist*(1/acc_score))*norm_dist
			##################
			spsp_paras.append((acc_score,diff_score,fgIC,bgIC,siteIC,dnadist,"".join(FGBP),"".join(BGBP)))
			spsp_sites.append(i+1)
			spsp_scores.append(sigvalue)
			usedsps=[ai for ai in splist if ai in sp_order]
	#get alignment information
	varsites=len(spsp_scores)
	return (gid,alignment_length,varsites,spsp_sites,spsp_paras,spsp_scores)

def batch_SpSpSite2(alidir,species,treefile,distfile,bgfile=None,fraction=1,gaps=1):
	alifilelist=glob("%s/*.fas"%alidir)
	alifilelist=[f for f in alifilelist if not f.endswith("anc.fas")]
	(DNAorder,distdat)=(parse_dist(distfile))
	sumspinfo=[]
	raw_dist, norm_dist, max_dist = calc_group_distance_normalized(treefile, species)
	for alifile in alifilelist:
		spspinfo=SpSpSite(alifile,species,norm_dist,DNAorder,distdat,bgfile=None,fraction=1,gaps=1)
		sumspinfo.append(spspinfo)
	return sumspinfo

def evo_acc(alidir,species,treefile,distfile,bgfile=None,STATFILE=sys.stdout,SUMFILE=sys.stdout,fraction=1,gaps=1):
	from scipy.stats import hypergeom 
	#print("gid\tsite\tJscore\tfgIC\tbgIC\tsiteIC\tdnadist\tICscore\tacc_score\tpvalue\tsignal\tbgAA\tfgAA",file=STATFILE)
	#print("alifile\tseqlen\tvar_sites\tsig_sites\tsum_score\tmax_score\tave_score\tpvalue\tacc",file=SUMFILE)
	sumspinfo=batch_SpSpSite2(alidir,species,treefile,distfile,bgfile=None,fraction=1,gaps=1)
	allvarsites=0
	allsigsites=[]
	allscores=[]
	nextsumpinfo=[]
	for (gid,alignment_length,varsites,spsp_sites,spsp_paras,spsp_scores) in sumspinfo:
		allscores+=spsp_scores+[0]*(alignment_length-varsites)
		allvarsites+=varsites
	data=np.array(allscores)
	"""
	from distfit import distfit
	dfit = distfit(distr=['norm','t','laplace','cauchy','chi2','expon','gamma','lognorm','uniform'])
	dfit.fit_transform(dacne)
	print("Best fit distribution:", dfit.model)
	"""
	if len(data)<1:
		print("%s contained no useable CNEs"%alidir)
		return 1
	score_mean=np.mean(data)
	shape, loc, scale = gamma.fit(data)
	fitted_gamma = gamma(a=shape, loc=loc, scale=scale)
	for (gid,alignment_length,varsites,spsp_sites,spsp_paras,spsp_scores) in sumspinfo:
		genepvalues=[]
		if len(spsp_scores)<1:
			continue
		for site,(acc_score,diff_score,fgIC,bgIC,siteIC,dnadist,FGBP,BGBP),score in zip(spsp_sites,spsp_paras,spsp_scores):
			pvalue=1 - fitted_gamma.cdf(score)
			genepvalues.append(pvalue)
			if pvalue < SIG_THRESHOLD:
				allsigsites.append(pvalue)
			print("{gid}\t{site}\t{J:.3f}\t{fgIC:.3f}\t{bgIC:.3f}\t{siteIC:.3f}\t\
{dnadist:.3f}\t{ICscore:.3f}\t{score:.6f}\t{pvalue:.6f}\t{signal}\t{bgbp}\t{fgbp}".format(
			gid=gid,
			site=site,
			J=diff_score,
			fgIC=fgIC,
			bgIC=bgIC,
			siteIC=siteIC,
			dnadist=dnadist,
			ICscore=acc_score,
			score=score,			
			pvalue=pvalue,
			signal = "ACC_POS" if score > score_mean and pvalue < SIG_THRESHOLD else \
				"ACC_NEG" if score < score_mean and pvalue > (1-SIG_THRESHOLD) else \
				"CON",
			bgbp=BGBP,
			fgbp=FGBP
			),file=STATFILE)
		nextsumpinfo.append((gid,alignment_length,varsites,spsp_scores,genepvalues))
		
	allsig=len(allsigsites)
	for (gid,alignment_length,varsites,spsp_scores,genepvalues) in nextsumpinfo:
		#output mafffiles
		#print(gid,alignment_length,spsp_scores)
		spsp_score_sum=sum(spsp_scores)
		spsp_score_max=max(spsp_scores)
		spsp_score_ave=spsp_score_sum/varsites
		#ave_pvalue=1 - fitted_gamma.cdf(spsp_score_ave)
		#p_gene = 1 - binom.cdf(k-1, m, alpha)
		sigpvalues=[ p for p in genepvalues if p < SIG_THRESHOLD ]
		n_sig=len(sigpvalues)
		#p_gene = 1 - binom.cdf(n_sig-1, varsites, SIG_THRESHOLD)
		#p_gene =1 - hypergeom.cdf(k-1, M, K, m)
		p_gene =1 - hypergeom.cdf(n_sig-1, allvarsites, allsig, varsites) 
		print("{gid}\t{seqlen}\t{varsites}\t{sigsites}\t{sum_score:.6f}\t{max_score:.6f}\t{ave_score:.6f}\t{pvalue:.6f}\t{acc}".format(
			gid=gid,
			seqlen=alignment_length,
			varsites=varsites,
			sigsites=n_sig,
			sum_score=spsp_score_sum,
			max_score=spsp_score_max,
			ave_score=round(spsp_score_ave,4),
			pvalue=p_gene,
			acc="ACC_POS" if spsp_score_ave > score_mean and p_gene < SIG_THRESHOLD else \
				"ACC_NEG" if spsp_score_ave < score_mean and p_gene > (1-SIG_THRESHOLD) else \
				"CON"
			),file=SUMFILE)
		
if __name__ == "__main__":
	alidir="bed_fasta/cne_500_200_diffr_CaseA"
	species="Tele,Pgut,Tbai,Nscu,Ptex,Nnaj"
	distfile="./dna.dist"
	treefile="./reptile.tre"
	SUMFILE=open("spsum.txt","w")
	STATFILE=open("statsum.txt","w")
	print("gid\tsite\tJscore\tfgIC\tbgIC\tsiteIC\tdnadist\tICscore\tacc_score\tpvalue\tsignal\tbgAA\tfgAA",file=STATFILE)
	print("alifile\tseqlen\tvar_sites\tsig_sites\tsum_score\tmax_score\tave_score\tpvalue\tacc",file=SUMFILE)
	Evo_acc(alidir,species,treefile,distfile,bgfile=None,STATFILE=STATFILE,SUMFILE=SUMFILE,fraction=1,gaps=1)
	STATFILE.close()
	SUMFILE.close()
