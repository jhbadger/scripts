#!/usr/bin/env python

from Bio import SeqIO
import bz2
import gzip
import argparse

# based on Trifonov EN, 1990 "Making sense of the human genome"
# calculates seq complexity based on words of size 1 to n seen
# versus total number of words of size 1 to n
def CT(seq, n):
    ct = 1
    for i in range(n):
        m = i+1
        vmax = float(min(4**m, len(seq)-m+1))
        words = dict()
        for j in range(0, len(seq), m):
            word = seq[j:j+m]
            words[word] = True
        ct *= len(words)/vmax
    return(ct)

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", required=True,
                                        help="input FASTQ file")
parser.add_argument("-t", "--threshold", default=0.1, type=float,
                                        help="complexity threshold")
parser.add_argument("-p", "--printval", action='store_true',
                    help="print seq and value")
args = parser.parse_args()

if ".bz2" in args.input:
    handle = bz2.BZ2File(args.input)
elif ".gz" in args.input:
    handle = gzip.open(args.input)
else:
    handle = open(args.input)

if args.printval:
    print "seq\tct"
    
for seq in SeqIO.parse(handle, "fastq"):
    ct = CT(str(seq.seq), 3)
    if args.printval:
        print seq.seq+"\t"+str(ct)
    elif ct >= args.threshold:
        print seq.format("fastq").rstrip("\n\r")
handle.close()
