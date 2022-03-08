#!/usr/bin/env python
# encoding: utf-8
"""
bed_from_genbank.py

grab the gene records from a genbank file (edit for other record types).

- requires:  biopython

"""

from Bio import SeqIO
import sys

def gb2bed(infile):
    print("\t".join(["Gene_ID","Gene_Type","Contig","Start","Stop","Strand","Gene_Name","Gene_Product"]))
    for record in SeqIO.parse(open(infile, "r"), "genbank") :
        #print(record.__dict__)
        chrom = record.name
        for feature in record.features:
            #print(feature.__dict__)
            if feature.type == 'CDS' or feature.type == "tRNA" or feature.type == "rRNA":
                start = feature.location.start.position
                stop = feature.location.end.position
                name = feature.qualifiers['locus_tag'][0]
                if "gene" in feature.qualifiers.keys():
                    gene = feature.qualifiers['gene'][0]
                else:
                    gene = ""
                if "product" in feature.qualifiers.keys():
                    product = feature.qualifiers['product'][0]
                else:
                    product = ""
                if feature.strand < 0:
                    strand = "-"
                else:
                    strand = "+"
                print("\t".join([name, feature.type, chrom, str(start),
                                 str(stop), strand, gene, product]))


if __name__ == '__main__':

    infile = ''
    try:
        infile = sys.argv[1]
    except:
        sys.exit("gb2bed.py INFILE")
    gb2bed(infile)
