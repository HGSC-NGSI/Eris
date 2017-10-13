#!/usr/bin/env python3

# Convert fluidigm raw csv to birdseed
# Based on csv2birdseed_v3.pl

from collections import defaultdict
import re
import sys


# Base compliment
def base_comp(seq):
    return seq.translate({ \
        ord('A'): 'T', \
        ord('C'): 'G', \
        ord('G'): 'C', \
        ord('T'): 'A'  \
    })


# Filter criteria for samples:
# num passing probe sites < 0.89 * 96 ~= 85
def filtered(sample):
    return len(sample) < 85


def run(probelist_fn, raw_bs_fn, ru2rs_fn):
    probes_dict = {}
    ru2rs_dict = {}
    birdseeds_dict = defaultdict(list)

    # Exclude from sample name characters that aren't alphanumeric, whitespace, or - _
    sample_name_re = re.compile(r'[^\d+\w+\-\_\s+]')
    # Match valid genotypes
    valid_geno_re = re.compile(r'[ACGT]:[ACGT]')
    # Split fields in files on , \t
    split_re = re.compile(r'[,\t]')

    # probes_dict {rsid: (chr, pos, ref, var), ... }
    print('Reading {}'.format(probelist_fn), file=sys.stderr)
    with open(probelist_fn) as probelist_fp:
        for line in probelist_fp:
            line = line.strip()
            fields = line.split('\t')
            probes_dict[fields[2]] = (fields[0], fields[1], fields[5], fields[6])
    print('Probelist size {}'.format(len(probes_dict.keys())), file=sys.stderr)

    # ru2rs_dict: {ruid: rsid, ... }
    print('Reading {}'.format(ru2rs_fn), file=sys.stderr)
    with open(ru2rs_fn) as ru2rs_fp:
        for line in ru2rs_fp:
            line = line.strip()
            fields = split_re.split(line)
            ru2rs_dict[fields[1]] = fields[0]
    print('RUID2RSID size {}'.format(len(ru2rs_dict.keys())), file=sys.stderr)

    print('Reading raw array file {}'.format(raw_bs_fn), file=sys.stderr)
    with open(raw_bs_fn) as raw_bs_fp:
        nr = 0

        for line in raw_bs_fp:
            nr += 1
            # skip 16 line header
            if nr <= 16:
                continue

            line = line.strip()
            fields = split_re.split(line)

            # Match valid genotypes
            if not valid_geno_re.match(fields[9]):
                continue

            # Birdseed filename
            fields[4] = sample_name_re.sub('', fields[4])
            sample_name = fields[4].split(' ')[0]
            uniq = fields[0].split('-')[0]
            birdseed_fn = '{}_{}.birdseed'.format(sample_name, uniq)

            ruid = fields[1]
            allele1, allele2 = fields[9].split(':')

            # Convert ruid to rsid (stored in ruid)
            if ruid in ru2rs_dict:
                ruid = ru2rs_dict[ruid]

            # birdseeds_dict:
            # {birdseed_fn: [[chr, pos, ref, allele1, allele2], ... ], ...}
            if ruid in probes_dict:
                chr, pos, ref, var = probes_dict[ruid]
                # from csv2birdseed.pl
                if fields[1] in ['hu202', 'hu317', 'hu44', 'hu68', 'hu69'] or \
                        allele1 not in [ref, var] or \
                        allele2 not in [ref, var]:
                    allele1 = base_comp(allele1)
                    allele2 = base_comp(allele2)
                birdseeds_dict[birdseed_fn].append([chr, pos, ref, allele1, allele2])

    # Write to birdseed files
    for fn in birdseeds_dict:
        if not filtered(birdseeds_dict[fn]):
            with open(fn, 'w') as fp:
                for entry in birdseeds_dict[fn]:
                    print('{}\t{}\t{}\t{}{}'.format(*entry), file=fp)
        else:
            print('{} filtered, num passing sites = {}'.format( \
                        fn, len(birdseeds_dict[fn])), file=sys.stderr)


if __name__ == '__main__':
    if len(sys.argv) < 4:
        print('python csv2birdseed.py probelist.txt raw_birdseed.csv ru2rs.txt', \
                file=sys.stderr)
        sys.exit(1)

    probelist_fn = sys.argv[1]
    raw_bs_fn = sys.argv[2]
    ru2rs_fn = sys.argv[3]

    run(probelist_fn, raw_bs_fn, ru2rs_fn)
