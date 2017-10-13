import sys
import os

#collect the probelist so I can get reference bases...
lookup = {}
probelist = sys.argv[1]
with open(probelist, 'r') as fh:
    for line in fh.readlines():
        data = line.strip().split('\t')
        lookup[data[2]] = [data[5],data[6]]
"""
Assumes
[Header]
GSGT Version	1.9.4
Processing Date	9/4/2015 1:42 PM
Content		humanomniexpress-24v1-0_a.bpm
Num SNPs	716503
Total SNPs	716503
Num Samples	24
Total Samples	24
File 	17 of 24
[Data]
Sample ID	SNP Name	Chr	Position	Allele1 - Design	Allele1 - Top	Allele2 - Top	SNP	Allele1 - Plus	Allele2 - Plus	GC Score	B Allele Freq	Log R Ratio
"""

missing = 0
with open(sys.argv[2], 'r') as fh:
    line = fh.readline().strip()
    while line != "[Data]":
        line = fh.readline().strip()
    line = fh.readline() # Header
    line = fh.readline() # process the first to get sample name
    data = line.strip().split('\t')
    with open(data[0] + ".birdseed", 'w') as fout:
        while line != "":
            if data[1] in lookup:
                fout.write("{chrom}\t{pos}\t{ref}\t{a1}{a2}\n" \
                    .format(chrom=data[2], pos=data[3], ref=lookup[data[1]][0], \
                    a1=data[8], a2=data[9]))
            else:
                missing += 1
            line = fh.readline()
            data = line.strip().split('\t')
print("%d missing rsid" % (missing))
