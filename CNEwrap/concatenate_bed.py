#!/usr/bin/env python3
import sys, os, re, argparse
from tqdm import tqdm


def parse_args():
    parser = argparse.ArgumentParser(
        description="Concatenate multiple FASTA files and output BED intervals."
    )
    parser.add_argument("-a", "--allsamples", required=True,
                        help="File containing list of all sample IDs (one per line).")
    parser.add_argument("-f", "--fastalist", required=True,
                        help="File containing list of fasta files (one per line).")
    parser.add_argument("-o", "--outprefix", required=True,
                        help="Output file prefix (e.g. concatenated).")
    return parser.parse_args()

def fas_merge(fasdict1, fasdict2):
    global fullkeys
    allfasdict = {}
    # pad Ns if missing
    for actkey in set(fullkeys) - set(fasdict1.keys()):
        fasdict1[actkey] = 'N' * len(next(iter(fasdict1.values())))
    for actkey in set(fullkeys) - set(fasdict2.keys()):
        fasdict2[actkey] = 'N' * len(next(iter(fasdict2.values())))
    assert sorted(fasdict1.keys()) == sorted(fasdict2.keys()), \
        'Error#: titles unequal between files'
    for fastitle in fasdict1:
        allfasdict[fastitle] = fasdict1[fastitle] + fasdict2[fastitle]
    return allfasdict

def parse_fas(fasfile):
    inhandle = open(fasfile)
    #sys.stderr.write('Parsing {0}...\n'.format(fasfile))
    fasdict = {}
    fastitle = re.compile(r"^>([\w.-]+)(\||\s*)")
    while True:
        lin = inhandle.readline()
        if not lin:
            break
        if lin.strip() == "":
            continue
        if lin[0] == ">":
            break
    while True:
        title = fastitle.match(lin).group(1)
        lin = inhandle.readline()
        while True:
            if not lin:
                break
            if lin[0] == ">":
                break
            fasdict[title] = fasdict.setdefault(title, "") + \
                             lin.strip().replace(" ", "").replace("\r", "")
            lin = inhandle.readline()
        if not lin:
            return fasdict

def main():
    args = parse_args()

    # read all sample IDs
    with open(args.allsamples) as samplelist:
        global fullkeys
        fullkeys = [s.strip() for s in samplelist if s.strip()]

    # read fasta file list
    with open(args.fastaList if hasattr(args, "fastaList") else args.fastalist) as flist:
        filelist = [line.strip() for line in flist if line.strip()]

    # parse all fasta
    allfas = []
    for fasfile in tqdm(filelist, desc="Parsing fasta files", unit="file"):
        allfas.append(parse_fas(fasfile))
    #allfas = [parse_fas(fasfile) for fasfile in filelist]

    # build concatenated
    fasdict = allfas[0]
    for d in allfas[1:]:
        fasdict = fas_merge(fasdict, d)

    # output concatenated fasta
    fas_out = args.outprefix + ".fas"
    with open(fas_out, "w") as fout:
        for faskey in sorted(fasdict):
            fout.write(">%s\n%s\n" % (faskey, fasdict[faskey]))

    # output bed file
    bed_out = args.outprefix + ".bed"
    with open(bed_out, "w") as bedout:
        for faskey in sorted(fasdict):
            start = 0
            for i, fasfile in enumerate(filelist):
                seglen = len(allfas[i][faskey])
                end = start + seglen
                bedout.write("{0}\t{1}\t{2}\t{3}\n".format(
                    faskey, start, end, os.path.basename(fasfile)))
                start = end

if __name__ == "__main__":
    main()
