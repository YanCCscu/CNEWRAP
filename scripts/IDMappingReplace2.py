#!/usr/bin/env python
from __future__ import print_function
import sys,re
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-m", "--mapfile", help="id_mapping file", action = "store")
parser.add_argument("-a", "--aimfile", help="aim file to revised", action = "store")
parser.add_argument("-o", "--outfile", help="prefix of output file", action = "store")
parser.add_argument("-r", "--reverse", help="replace items from colume One with corresponding items in colume Two in mapfile",\
			action = "store_true", default = False)
args = parser.parse_args()
outfile=open(args.outfile,'w')
#restring=sys.argv[3]
with open(args.aimfile,'r') as targetfile:
	strline = "".join(targetfile.readlines())
	with open(args.mapfile,'r') as mapfile:
		if args.reverse:
			for l in mapfile:
				strline=strline.replace(l.split()[1],l.split()[0]) 
		else:
			for l in mapfile:
				strline=strline.replace(l.split()[0],l.split()[1])
		print(strline,end="",file=outfile)

outfile.close()
		
