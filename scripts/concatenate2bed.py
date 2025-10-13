#!/usr/bin/env python
from functools import reduce
from glob import glob
import sys,os,re
introduction="\n\
------------------------------------------------------------------\n\
This program align sequences in you files base on those equal names\n\
or the names with same begin:\n\
USAGE:python concatenate.py infile1 infile2 infile3 ... \n\
      python concatenate.py dir ... \n\
NOTE: the id should match in each files, the first strings after \'>\'\n\
and before spaces and \'|\' are assigned id of its following sequence\n\
--------------------------------------------------------------------\n"
if len(sys.argv)<2:
	print ('\033[1;32;40m',)	#make sure there are imput files.
	print (introduction)
	print ('\033[0m',)
	sys.exit(1)
"""---------------------------------------------------------------------	"""
def fas_merge(faslist):      #merge 2 fasta file contain same set of titles
	global fullkyes
	allfasdict={}
	seqids=[]
	seqlens=[]
	for fasfile in faslist:
		fasdict=parse_fas(fasfile)
		seqlen=len(next(iter(fasdict.values())))
		for actkey in set(fullkeys)-set(fasdict.keys()):
			fasdict[actkey]='N'*seqlen
		for fastitle in fasdict:
			allfasdict.setdefault(fastitle,[]).append(fasdict[fastitle])
		seqids.append(os.path.basename(fasfile).replace(".fas","").replace(".fasta","").replace(".fa",""))
		seqlens.append(seqlen)
	return (allfasdict,seqids,seqlens)
	
	
def parse_fas(fasfile): #parse fasta and return a dict 
	inhandle=open(fasfile)
	#sys.stderr.write('Parsing {0}...\n'.format(fasfile))
	fasdict={}
	fastitle=re.compile(r"^>([\w.-]+)(\||\s*)")
	while True:                #skip empty lines and get the first title line
		lin=inhandle.readline()
		if not lin:
			break
		if lin=="":
			continue
		if lin[0]==">":
			break
	while True:		#compose the dict storing title and sequence
		title=fastitle.match(lin).group(1)
		lin=inhandle.readline()
		while True:
			if not lin:
				break
			if lin[0]==">":
				break
			fasdict[title]=fasdict.setdefault(title,"")+lin.strip().replace(" ","").replace("\r","")
			lin=inhandle.readline()
		if not lin:
			#sys.stderr.write('Finished parsing %s\n'%fasfile)
			return fasdict
def cumul(lst):
    product = 0
    result = []
    for num in lst:
        product += num
        result.append(product)
    return result
"""-------------------------------------------------------------------------------"""
sampefile=sys.argv[1]
with open(sampefile) as samplelist:
	fullkeys=[s.strip() for s in samplelist]
if os.path.isdir(sys.argv[2]):
	alilist=glob("%s/*.fas"%sys.argv[2])
	filelist=[ bedfas for bedfas in alilist if not ".anc.fas" in bedfas]
	keywd=os.path.basename(sys.argv[2])
else:
	filelist=sys.argv[2:]
	keywd="All"
(fasdict,ids,lens)=fas_merge(filelist)
#print(fasdict)
#print(iddict)
with open(keywd+"_Cated_seqs.fas",'w') as FA:
	for faskey in sorted(fasdict):
		print (">%s\n%s"%(faskey,"".join(fasdict[faskey])),file=FA)
with open(keywd+"_Cated_seqs.bed",'w') as BED:
	starts=cumul([0]+lens[:-1])
	ends=cumul(lens)
	for seqid,start,end in zip(ids,starts,ends):
		print("{chr}\t{start}\t{end}\t{id}".format(chr=fullkeys[0],start=start,end=end,id=seqid),file=BED)
