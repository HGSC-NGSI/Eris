#!/usr/bin/python

"""Eris is a program for calculating sample concordance between BAM/CRAM
alignment files and SNP array files in birdseed format. Also provided are
scripts for converting array files for some common SNP chips (such as
Fluidigm SNPTrace and OmniExpress1.0) into the birdseed format.

Contributors:
    Adam English, Adam Mansfield, Shruthi Ambreth, Jeffrey Reid

Copyright 2017 Baylor College of Medicine Human Genome Sequencing Center
"""

import argparse
import glob
import logging
import math
import os
import pickle
import pysam
import re
import sys
from collections import defaultdict, OrderedDict


USAGE = r"""
Estimate concordance and contamination of BAMs/CRAMs against SNP array 
files in birdseed format.

Tested on pysam 8.4 and python 3.5.

Running Eris for one BAM:
    > python3 Eris.py -b BAM -a ARRAY_DIR -A SELF_ARRAY -p PROBELIST -o OUTFILE
    
    e.g. for Fluidigm SNPTrace arrays:
    > python3 Eris.py -b NA12878.bam -a /path/to/birdseeds_dir \
      -A /path/to/birdseeds_dir/NA12878.birdseed -o NA12878.report \
      -p ../array_kits/fluidigm_37/probelist.txt
    
Running Eris for many BAMs (mostly specific to HGSC):
    If the input is --user_input, we expect a text file with the following 
    tab-delimited columns (no header row):
        run_name   - identifier for the job (can be none)
        sample     - Sample identifier (must match the SM in BAM's header)
        result_path - Path to a directory (Sequencing or Merge) with a BAM
        array_type  - type of array
        array_dir   - Path to raw array data
    
    This will return a bash script with the commands to run an Eris job for
    each line of input on the MOAB cluster.
"""


def setup_logging(debug=False):
    log_level = logging.DEBUG if debug else logging.INFO
    log_format = "%(asctime)s [%(levelname)s] %(message)s"
    logging.basicConfig(stream=sys.stderr, level=log_level, format=log_format)
    logging.info("Running %s" % " ".join(sys.argv))


def parse_args():
    parser = argparse.ArgumentParser(prog="Eris", description=USAGE,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-u", "--user_input", type=str, default=None,
                        help="Create commands from user_input for execution (Ignores)")
    parser.add_argument("-b", "--bam", metavar="BAM/CRAM", type=str,
                        help="Sample's BAM or CRAM file. Either must be indexed.")
    parser.add_argument("-p", "--probelist", metavar="PROBE", type=str,
                        help="Probe definition file (1-based coordinates)")
    parser.add_argument("-a", "--array_dir", metavar="ARRAYS", type=str,
                        help="Directory with all arrays to compare against.")
    parser.add_argument("-A", "--Array", type=str,
                        help="Force this as sample's self array")
    parser.add_argument("-o", "--output", default=None,
                        help="Name of output file (stdout)")
    parser.add_argument("--debug", action="store_true",
                        help="Verbose logging")
    parser.add_argument("-g", "--genotype_output",
                        help="Preprocess a BAM into a file containing only "
                             "genotypes that Eris needs to run concordance "
                             "with a given probelist. If you aren't sure "
                             "which probelist you will need in the future, "
                             "give a probelist containing all unique sites "
                             "from probelists you might need. Provided arg "
                             "to -g is location of output.")
    parser.add_argument("-i", "--in_silico",
                        help="Path to a pkl genotype file created from a bam "
                             "using --genotype_output. This can be used instead "
                             "of an input bam.")
    parser.add_argument("-f", "--force", action="store_true",
                        help="Creates commands from user_input even if no "
                             "matching array is present.")
    parser.add_argument("-F", "--fake_arrays", action="store_true",
                        help="Creates commands from user_input even if no "
                             "array is present and creates fake arrays "
                             "using the sample's name.")

    args = parser.parse_args()
    setup_logging(args.debug)
    if args.bam is None and args.user_input is None and args.in_silico is None:
        logging.error("Must provide --bam, --user_input, or --in_silico")
        print("Must provide --bam, --user_input, or --in_silico")
        exit(1)

    if args.array_dir is None and args.genotype_output is None:
        logging.error("Must provide --array_dir or --genotype_output")

    if args.output is None:
        args.output = sys.stdout
    else:
        args.output = open(args.output, 'w')

    return args


# NOTE: Coordinates in pysam are always 0-based (following the python
# convention). SAM/BAMs use 1-based coordinates.

# Arbitrarily defined by users
MAX_CUTOFF = 10000  # maximum number of reads to call genotype
MIN_CUTOFF = 4  # minimum number of reads to call genotype
HOM_REF_CUTOFF = 0.1  # alt_reads < cutoff * total_reads means homozygous ref
HOM_ALT_CUTOFF = 0.75  # alt_reads > cutoff * total_reads means homozygous alt

########################################
# -- Nucleotide to Genotype Helpers -- #
########################################

NUCS = {"A": 1,
        "T": 2,
        "C": 4,
        "G": 8,
        "N": 16}

GTS = dict((y, x) for x, y in NUCS.items())


def nuc2gt(nuc):
    """Turn allele(s) into genotype NUCS code."""
    a, b = list(nuc)
    if a == b:
        return NUCS[a]
    return NUCS[a] + NUCS[b]


def gt2nuc(gt):
    """Turn genotype NUCS code into list of nucleotides."""
    if gt in GTS:
        return [GTS[gt], GTS[gt]]
    return [GTS[x] for x in GTS if gt & x]


def gt_compare(self, other):
    """Number of alleles common between self and other."""
    if self == other:
        return 2
    elif self.gt & other.gt:
        return 1

    return 0


def is_gt_homozygous(gt):
    """
    Return whether genotype is homozygous.

    Args:
        gt (int):  genotype NUCS code
    Returns:
        bool
    """
    return gt in GTS


class ArraySite(object):
    """Represents one site in an array."""

    def __init__(self, ref, genotype, nul_count=0, ref_count=0, alt_count=0):
        self.ref = ref
        self.genotype = genotype
        self.null_count = nul_count
        self.ref_count = ref_count
        self.alt_count = alt_count

    def to_tuple(self):
        return self.ref, self.genotype, self.null_count, self.ref_count, self.alt_count


class ProbeSite(object):
    """Represents a probelist site for an array.
    Will also hold sample ref/alt/genotype for whatever we're holding.
    """

    def __init__(self, chrom, pos, rsid, left_flank, right_flank, ref_allele, alt_allele,
                 ref_freq, het_freq, var_freq):
        """
        Args:
            chrom (str)
            pos (int)
            rsid (str)
            left_flank (str)
            right_flank (str)
            ref_allele (str)
            alt_allele (str)
            ref_freq (float)
            het_freq (float)
            var_freq (float)
        """
        self.chrom = chrom
        self.pos = int(pos)
        self.rsid = rsid
        self.left_flank = left_flank
        self.right_flank = right_flank
        self.ref_allele = ref_allele
        self.alt_allele = alt_allele
        self.ref_freq = float(ref_freq)
        self.het_freq = float(het_freq)
        self.var_freq = float(var_freq)

    def __repr__(self):
        values = ", ".join([str(x) for x in [
            self.left_flank,
            self.right_flank,
            self.ref_allele,
            self.alt_allele,
            self.ref_freq,
            self.het_freq,
            self.var_freq
            ]])
        return "<Probesite %s at (%s:%d): %s>" % (self.rsid, self.chrom, self.pos, values)

    def __str__(self):
        return ", ".join([str(x) for x in [
            self.chrom,
            self.pos,
            self.rsid,
            self.left_flank,
            self.right_flank,
            self.ref_allele,
            self.alt_allele,
            self.ref_freq,
            self.het_freq,
            self.var_freq
            ]])


class ProbeBamReport(object):
    def __init__(self):
        self.num_probe_sites = 0  # How many probe sites we were given and had valid reference to lookup
        self.num_bad_cov_sites = 0  # How many sites we couldn't genotype
        self.num_homref_sites = 0  # How many reference homozygous sites we're filtering out
        self.num_het_sites = 0  # How many sites we genotyped as het
        self.num_homalt_sites = 0  # How many sites we genotyped as hom alt
        self.total_reads_analyzed = 0  # Total number of reads we've looked at

    def pretty_report(self):
        s = "##ProbeSitesCount        %d\n" % self.num_probe_sites
        s += "##CoverageFiltStes       %d\n" % self.num_bad_cov_sites
        s += "##HomRefSitesCount       %d\n" % self.num_homref_sites
        s += "##HetSitesCount          %d\n" % self.num_het_sites
        s += "##HomAltSitesCount       %d\n" % self.num_homalt_sites
        s += "##NumReadsAnalyzed       %d" % self.total_reads_analyzed
        return s


class ContaminationReport(object):
    def __init__(self):
        self.sites = 0  # Number of sites where we can look for contamination
        self.site_reads = 0  # Number of reads at contamination sites
        self.null_reads = 0  # Reads that are contaminated
        self.match_reads = 0  # Reads that match the expected homozygous allele
        self.miss_reads = 0  # Reads that are 'contamination' but actually match the non-homozygous allele
        self.null_score = 0  # Reads that are contaminated
        self.match_score = 0  # Reads that match the expected homozygous allele
        self.miss_score = 0  # Reads that are 'contamination' but actually match the non-homozygous allele

        self.contamination = -1

    def calc_contamination(self):
        if self.site_reads == 0:
            self.contamination = 0.0
        else:
            numerator = self.null_score + self.miss_score
            denominator = numerator + self.match_score
            # self.contamination = (float(self.null_reads + self.miss_reads) / self.site_reads) * 100
            self.contamination = (numerator / denominator) * 100

    def pretty_report(self):
        s = "##ContaminationSites     %d\n" % self.sites
        s += "##ContaminationSiteReads %d\n" % self.site_reads
        s += "##ContaminationMatReads  %d\n" % self.match_reads
        s += "##ContaminationMisReads  %d\n" % self.miss_reads
        s += "##ContaminationNulReads  %d\n" % self.null_reads
        s += "##Contamination          %.2f%%" % self.contamination
        return s


class ConcordanceReport(object):
    def __init__(self):
        self.array_sites = 0  # number of sites in the array
        self.array_hom = 0  # array homozygous sites
        self.matched_sites = 0  # number of sites matched between in_silico and array
        self.unmatched_sites = 0  # sites we couldn't match up
        self.concordant_num = 0  # concordant meaning at least one of the alleles match
        self.nonconcordant_num = 0  # non-concordant meaning none of the alleles match

        self.exact_match_count = 0  # sequence matches array
        self.exact_match_hom = 0  # score of homozygous sites
        self.exact_match_het = 0  # score of heterozygous sites

        self.one_match_count = 0  # part of sequence matches array
        self.one_match_score = 0
        self.one_mismatch_score = 0

        self.no_match_score = 0

        self.concordance = 0

    def calc_concordance(self):
        """This is the concordance formula -- essentially concordant alleles
        over all alleles checked.

        Returns:
            float
        """
        numerator = self.exact_match_het * 2 + self.exact_match_hom * 2 + self.one_match_score
        denominator = self.one_match_score + self.one_mismatch_score + self.exact_match_het * 2 + \
            self.exact_match_hom * 2 + self.no_match_score * 2

        if denominator == 0:
            self.concordance = 0.0
        else:
            self.concordance = (float(numerator) / float(denominator)) * 100
        return self.concordance

    @staticmethod
    def get_header():
        """Get \t delimited string of the header"""
        return "\t".join([
            "ArraySites",
            "ArrayHomSites",
            "MatchedSites",
            "UnmatchedSites",
            "ConcordantSites",
            "NonConcordantSites",
            "ExactMatchCnt",
            "ExactMatchHomScore",
            "ExactMatchHetScore",
            "OneMatchCount",
            "OneMatchScore",
            "OneMismatchScore",
            "NoMatchScore",
            "Concordance"
            ])

    def pretty_report(self):
        s = "##ArraySites             %d\n" % self.array_sites
        s += "##ArrayHomSites          %d\n" % self.array_hom
        s += "##MatchedSites           %d\n" % self.matched_sites
        s += "##UnmatchedSites         %d\n" % self.unmatched_sites
        s += "##ConcordantSites        %d\n" % self.concordant_num
        s += "##NonConcordantSites     %d\n" % self.nonconcordant_num
        s += "##ExactMatchCnt          %d\n" % self.exact_match_count
        s += "##ExactMatchHomScore     %.2f\n" % self.exact_match_hom
        s += "##ExactMatchHetScore     %.2f\n" % self.exact_match_het
        s += "##OneMatchCount          %d\n" % self.one_match_count
        s += "##OneMatchScore          %.2f\n" % self.one_match_score
        s += "##OneMisMatchScore       %.2f\n" % self.one_mismatch_score
        s += "##NoMatchScore           %.2f\n" % self.no_match_score
        s += "##Concordance            %.2f%%" % self.concordance

        return s

    def __str__(self):
        s = "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2f\t%.2f\t%d\t%.2f\t%.2f\t%.2f\t%.2f%%"
        return s % (self.array_sites,
                    self.array_hom,
                    self.matched_sites,
                    self.unmatched_sites,
                    self.concordant_num,
                    self.nonconcordant_num,
                    self.exact_match_count,
                    self.exact_match_hom,
                    self.exact_match_het,
                    self.one_match_count,
                    self.one_match_score,
                    self.one_mismatch_score,
                    self.no_match_score,
                    self.concordance)


def build_probe_sites(probe_definition):
    """Given a probe definition and array data, build all the ProbeSites as
    {chrom: {pos: ProbeSite}}

    Args:
        probe_definition (str): Path to probelist file

    Returns:
        dict[str, dict[int, ProbeSite]]
    """
    probe_count = 0
    ret = defaultdict(dict)

    # Use all probesites
    with open(probe_definition, 'r') as fh:
        for line in fh:
            probe = ProbeSite(*line.strip().split('\t'))
            ret[probe.chrom][probe.pos] = probe
            probe_count += 1
    logging.info("%d probe sites" % probe_count)
    return ret


def build_array(probe_sites, birdseed_filename):
    """
    Create array {chrom: {pos: ArraySite}} for birdseed_filename at probe_sites

    Args:
        probe_sites (dict[str, dict[int, ProbeSite]])
        birdseed_filename (str)

    Returns:
        dict[str, dict[int, ArraySite]]

    """
    ret = defaultdict(dict)

    with open(birdseed_filename, 'r') as fh:
        logging.debug(birdseed_filename)
        filter_count = 0

        if 'chr' in ''.join(list(probe_sites.keys())):
            filter_chr = False
        else:
            filter_chr = True

        for line in fh:
            try:
                chrom, pos, ref, alleles = line.strip().split()

                if filter_chr:
                    chrom = chrom.replace('chr', '')

                if '-' in alleles or '0' in alleles:
                    filter_count += 1
                    continue

                pos = int(pos)

                if chrom in probe_sites and pos in probe_sites[chrom]:
                    ret[chrom][pos] = ArraySite(ref, nuc2gt(alleles))

            except (ValueError, KeyError, AttributeError):
                logging.warning("Skipping weird site '{}' in {}".format(line.strip(), birdseed_filename))

    if filter_count > 0:
        logging.warning("Filtered {} sites in birdseed {} with bad alleles.".format(filter_count, birdseed_filename))

    return dict(ret)


def probe_bam(bam, probe_sites):
    """Given a bam and probe_sites, goes to each site in the bam and calls a
    genotype at that location.

    Args:
        bam (pysam.AlignmentFile):  bam to probe
        probe_sites (dict(str, dict(int, ProbeSite))):  sites to probe bam at

    Returns:
        dict[str, dict[int, ArraySite]:   dict of genotypes at each site
    """
    logging.info("Genotyping bam {} at probe sites".format(bam.filename))

    genotypes = defaultdict(dict)
    report = ProbeBamReport()

    for chrom in probe_sites:
        if chrom not in bam.references:
            logging.warning("Chromosome {} not in reference".format(chrom))
            continue

        for pos in probe_sites[chrom]:
            probe = probe_sites[chrom][pos]
            report.num_probe_sites += 1
            pos = probe.pos - 1  # SAM is 1-based, pysam is 0 based

            ref_count = 0
            alt_count = 0
            null_count = 0

            for pileup_col in bam.pileup(probe.chrom, pos, pos + 1):
                if pos != pileup_col.pos:
                    continue

                for pileup_read in pileup_col.pileups:
                    if pileup_read.indel or pileup_read.query_position is None:
                        continue

                    report.total_reads_analyzed += 1

                    base = pileup_read.alignment.seq[pileup_read.query_position]

                    if base == probe.ref_allele:
                        ref_count += 1
                    elif base == probe.alt_allele:
                        alt_count += 1
                    else:
                        null_count += 1

            total_bases = ref_count + alt_count

            if total_bases <= MIN_CUTOFF:
                logging.debug("Low coverage at {}:{}, no genotype possible".format(probe.chrom, probe.pos))
                genotype = None
            elif total_bases >= MAX_CUTOFF:
                logging.debug("High coverage at {}:{}, no genotype possible".format(probe.chrom, probe.pos))
                genotype = None
            elif alt_count < HOM_REF_CUTOFF * total_bases:  # hom ref
                logging.debug("Probe {}:{} is 0/0 {}".format(probe.chrom, probe.pos, probe.ref_allele))
                genotype = NUCS[probe.ref_allele]
                report.num_homref_sites += 1
            elif alt_count > HOM_ALT_CUTOFF * total_bases:  # hom alt
                logging.debug("Probe {}:{} is 1/1 {}".format(probe.chrom, probe.pos, probe.ref_allele))
                genotype = NUCS[probe.alt_allele]
                report.num_homalt_sites += 1
            else:  # het
                logging.debug("Probe {}:{} is 0/1 {}".format(probe.chrom, probe.pos, probe.ref_allele))
                genotype = NUCS[probe.ref_allele] + NUCS[probe.alt_allele]
                report.num_het_sites += 1

            if genotype is None:
                report.num_bad_cov_sites += 1
                continue

            genotypes[probe.chrom][probe.pos] = ArraySite(ref=probe.ref_allele, genotype=genotype, nul_count=null_count,
                                                          ref_count=ref_count, alt_count=alt_count)

    return dict(genotypes), report


def save_in_silico_array(in_silico_array, probe_bam_report, sample_name, filename):
    """Converts in-silico array from dict of dicts of ArraySites into a
    dict of dicts of tuples, then pickles it along with its associated
    ProbeBamReport at filename.

    Args:
        in_silico_array (dict[chrom, dict[pos, ArraySite]])
        probe_bam_report (ProbeBamReport)
        sample_name (str)
        filename (str)
    """
    tuple_in_silico_array = defaultdict(dict)

    for chrom in in_silico_array:
        for pos in in_silico_array[chrom]:
            tuple_in_silico_array[chrom][pos] = in_silico_array[chrom][pos].to_tuple()

    tuple_in_silico_array = dict(tuple_in_silico_array)

    with open(filename, 'wb') as pickle_file:
        pickle.dump((tuple_in_silico_array, probe_bam_report, sample_name), pickle_file)


def load_in_silico_array(filename):
    """Retrieve and format in-silico genotype array and associated
    ProbeBamReport and sample name stored in pickle file at given filename.

    Args:
        filename (str): Path to pickle file containing
            tuple(in_silico_genotypes, ProbeBamReport, sample_name)

    Returns:
        tuple(dict[str, dict[int, ArraySite]], ProbeBamReport, str)
    """
    with open(filename, 'rb') as pickle_file:
        tuple_in_silico_array, probe_bam_report, sample_name = pickle.load(pickle_file)

    in_silico_array = defaultdict(dict)

    for chrom in tuple_in_silico_array:
        for pos in tuple_in_silico_array[chrom]:
            in_silico_array[chrom][pos] = ArraySite(*tuple_in_silico_array[chrom][pos])

    in_silico_array = dict(in_silico_array)

    return in_silico_array, probe_bam_report, sample_name


def generate_contamination_report(in_silico_array, self_array, probe_sites):
    """Returns a Contamination report for an array genotyped from a BAM against
    its matching SNP array. Args specify {chrom: {pos: site}

    Args:
        in_silico_array (dict[str, dict[int, ArraySite]]):    Genotypes from BAM
        self_array (dict[str, dict[int, ArraySite]]):   Corresponding SNP array
        probe_sites (dict[str, dict[int, ProbeSite]]):  Sites to examine

    Returns:
        ContaminationReport
    """
    report = ContaminationReport()

    for chrom in probe_sites:
        for pos in probe_sites[chrom]:
            probe = probe_sites[chrom][pos]
            silico_site = in_silico_array.get(probe.chrom, {}).get(probe.pos, None)
            self_site = self_array.get(probe.chrom, {}).get(probe.pos, None)

            if silico_site is None or self_site is None:
                continue

            if self_site.genotype in GTS:  # then it's homozygous something
                report.sites += 1
                report.site_reads += silico_site.ref_count + silico_site.alt_count

                report.null_reads += silico_site.null_count
                report.null_score += silico_site.null_count * 0.01  # Error rate

                if GTS[self_site.genotype] == probe.ref_allele:
                    report.match_reads += silico_site.ref_count
                    report.match_score += silico_site.ref_count * (1 - probe.ref_freq)
                    report.miss_reads += silico_site.alt_count
                    report.miss_score += silico_site.alt_count * (1 - probe.var_freq)

    report.calc_contamination()

    return report


def site_compare(probe, site1, site2, report):
    """Takes a probe site and compares the two values from an in-silico array
    and a SNP array and adds to the Concordance report.

    site1 is assumed to be the in-silico array.

    Args:
        probe (ProbeSite)
        site1 (ArraySite):  Usually the site from the in-silico (BAM) array
        site2 (ArraySite):  Usually the site from the SNP array chip
        report (ConcordanceReport)
    """
    report.matched_sites += 1
    if site2.genotype in GTS:
        report.array_hom += 1
        if site2.genotype == NUCS[probe.ref_allele]:
            return
    # Exact match
    if site1.genotype == site2.genotype:
        report.concordant_num += 1
        report.exact_match_count += 1
        if site1.genotype in GTS:  # homozygous
            report.exact_match_hom += 1 - probe.var_freq
        else:  # het
            report.exact_match_het += 1 - probe.het_freq
    # Single match
    elif site1.genotype & site2.genotype:
        report.concordant_num += 1
        report.one_match_count += 1
        if site1.genotype in GTS:  # homozygous
            report.one_match_score += 1 - probe.var_freq
            report.one_mismatch_score += 1 - probe.var_freq
        else:  # het
            report.one_match_score += 1 - probe.het_freq
            report.one_mismatch_score += 1 - probe.het_freq
    else:  # No match
        report.nonconcordant_num += 1
        if site1.genotype in GTS:  # homozygous
            report.no_match_score += 1 - probe.var_freq
        else:
            report.no_match_score += 1 - probe.het_freq


def make_concordance(probe_sites, in_silico, sample_array):
    """Given a sample's sequencing bam probes and in_silico array, compare it
    against sample_array values and create a concordance report.

    Args:
        probe_sites (dict[str, dict[int, ProbeSite]])
        in_silico (dict[str, dict[int, ArraySite]])
        sample_array (dict[str, dict[int, ArraySite]])

    Returns:
        ConcordanceReport
    """
    con = ConcordanceReport()
    for chrom in probe_sites:
        if chrom not in sample_array:
            con.unmatched_sites += len(probe_sites[chrom])
            continue

        con.array_sites += len(sample_array[chrom])
        if chrom not in in_silico:
            con.unmatched_sites += len(sample_array[chrom])
            continue

        for pos in probe_sites[chrom]:
            if pos not in sample_array[chrom] or pos not in in_silico[chrom]:
                con.unmatched_sites += 1
                continue

            probe = probe_sites[chrom][pos]
            site1 = in_silico[chrom][pos]
            site2 = sample_array[chrom][pos]
            site_compare(probe, site1, site2, con)

    con.calc_concordance()

    return con


def output_reports(ret):
    """Sort concordance results based on concordance and calculate the average."""
    values = []
    total = 0.0
    for i in ret:
        conc = ret[i].concordance
        values.append((conc, i))
        total += conc

    values.sort(reverse=True)
    ret2 = OrderedDict()
    for conc, key in values:
        ret2[key] = ret[key]

    return ret2, total / len(ret2)


def batch_stats(reports, self_hit):
    """
    Args:
        reports (iterable[str, ConcordanceReport])
        self_hit (ConcordanceReport)

    Returns:
        (float, float): Mean and z-score
    """
    mean = sum(report.concordance for name, report in reports) / len(reports)
    if len(reports) - 1 == 0 or self_hit is None:
        return mean, 0

    sd = math.sqrt(sum([(mean - report.concordance) ** 2 for name, report in reports]) / (len(reports) - 1))
    if sd == 0:
        return mean, 0

    return mean, (self_hit.concordance - mean) / sd


def run_judgement(sorted_reports, self_hit, contamination):
    """Given sorted concordance reports, average concordance, and contamination,
    make a judgement.

    Args:
        sorted_reports (list[str, ConcordanceReport])
        self_hit (ConcordanceReport)
        contamination (ContaminationReport)

    Returns:
        str:    Judgement, e.g. "POSSIBLE CONTAMINATION"
    """
    avg, z = batch_stats(sorted_reports, self_hit)
    if len(sorted_reports) == 1:
        return "SINGLE ARRAY", avg, z

    if self_hit is None:
        return "NOSAMPLEARRAY", avg, z

    best_hit, second_hit = [x[1] for x in sorted_reports[:2]]

    if self_hit.concordant_num + self_hit.nonconcordant_num <= 500:
        # (but there is a fluidigm)
        return "INSENSITIVE", avg, z

    if avg < 50 or avg > 75:
        return "INCONCLUSIVE", avg, z

    contam = contamination.contamination

    if self_hit.concordance > 0.90:
        if contam > 0.05:
            return "POSSIBLE CONTAMINATION", avg, z

        if best_hit.concordance >= 0.90 and second_hit.concordance >= 0.90:
            return "POSSIBLE SAME SUBJECT MATCH", avg, z

        if self_hit == best_hit:
            return "PASS", avg, z

        return "POSSIBLE SAME SUBJECT MATCH", avg, z

    if self_hit.concordance > 0.80:
        if self_hit == best_hit:
            if contam > 0.05:
                return "POSSIBLE CONTAMINATION", avg, z
            else:
                return "MARGINAL CONCORDANCE", avg, z
        if best_hit.concordance <= 0.90:
            if contam > 0.05:
                return "POSSIBLE CONTAMINATION OR SWAP", avg, z
            else:
                return "POSSIBLE SAME SUBJECT MATCH", avg, z
        return "KNOWN SWAP OR SAME SUBJECT MATCH", avg, z

    if self_hit == best_hit:
        if contam > 5:
            return "POSSIBLE CONTAMINATION", avg, z
        return "UNDETERMINED", avg, z

    if best_hit.concordance <= 0.90:
        if contam > 0.05:
            return "POSSIBLE CONTAMINATION OR SWAP", avg, z
        return "UNDETERMINED", avg, z

    return "NO JUDGEMENT MADE", avg, z


def make_commands(args):
    """Just print the commands that we would create for every user_input line."""
    cmds = []
    with open(args.user_input, 'r') as fh:
        for line in fh:
            if '\t' in line.strip():
                run_name, sample, result_path, array_type, raw_array = line.strip().split('\t')
            else:
                run_name, sample, result_path, array_type, raw_array = line.strip().split()

            logging.info("Creating command for %s" % sample)

            out = sample + ".report"

            bam = find_bam(result_path)
            if bam is None:
                logging.error("No Bam Found for %s" % sample)
                continue

            ret = find_birdseed(sample, args.array_dir, args.Array)
            if ret is not None:
                my_array, all_arrays = ret
            else:
                my_array, all_arrays = (None, None)

            if my_array is None:
                if args.force:
                    pass
                elif args.fake_arrays:
                    fake_path = "{}_FAKE.birdseed".format(sample)
                    logging.warning("Creating fake birdseed {}".format(fake_path))
                    my_array = os.path.join(args.array_dir, fake_path)
                else:
                    continue

            cmd = "{pg} -b {bam} -p {probe} -a {array_dir} -A {my_array} -o {out}" \
                .format(pg=__file__, bam=bam, probe=args.probelist,
                        array_dir=args.array_dir, my_array=my_array, out=out)
            cmd = 'echo "{cmd}" | msub -N {sample}_eris -d `pwd`/ -l nodes=1:ppn=1,mem=8000mb' \
                .format(cmd=cmd, sample=sample)
            cmds.append(cmd)
    print("\n".join(cmds))


def find_bam(result_path):
    """Returns the single bam in the result path.

    Args:
        result_path (str):  Path to sequencing event (specific to HgV/HGSC).
                            Basically just looking in directory for bams or
                            crams.
    Returns:
        str or None
    """
    bams = glob.glob(os.path.join(result_path, "*.bam"))
    bams.extend(glob.glob(os.path.join(result_path, "*.cram")))
    if len(bams) != 1:
        logging.warning("Exactly one bam not found in %s" % result_path)
    for i in bams:
        if i.endswith(".hgv.bam"):
            return i
    for i in bams:
        if i.endswith(".bam"):
            return i
    for i in bams:
        if i.endswith(".hgv.cram"):
            return i
    for i in bams:
        if i.endswith(".cram"):
            return i
    return None


def find_birdseed(sample_name, array_dir, forced_array=None):
    """Returns the single array that matches or was forced in a tuple with all
    of the other arrays.

    Args:
        sample_name (str)
        array_dir (str)
        forced_array (dict[str, dict[int, ArraySite]] or None)

    Returns:
        ArraySite, list[ArraySite]
    """
    # need to find my raw array
    all_arrays = glob.glob(os.path.join(array_dir, "*.birdseed"))
    logging.debug(all_arrays)

    if forced_array is not None:
        if all_arrays.index(forced_array):
            all_arrays.remove(forced_array)
        return forced_array, all_arrays

    my_array = [x for x in all_arrays if re.search(r"\b%s\b" % sample_name.replace('-', '\-'), x) is not None]
    # Pin down which one is exactly my_array
    if len(my_array) != 1:
        logging.debug("No regex match")
        my_array = [x for x in all_arrays if x.count(sample_name)]
        logging.debug(my_array)
        if len(my_array) != 1:
            logging.warning("Couldn't find exactly one array matching sample name %s" % sample_name)
            logging.warning(my_array)

            if my_array:
                # Coose youngest array, prefer arrays that aren't conversions (idiomatic to HGSC)
                no_conv = [m for m in my_array if 'conv' not in m]
                if no_conv:
                    my_array = no_conv
                my_array = max(my_array, key=os.path.getmtime)
                logging.warning("Choosing %s by default" % my_array)
            else:
                my_array = None
        else:
            if 'FAKE' in my_array[0]:
                logging.warning("Couldn't find exactly one array matching sample name %s" % sample_name)
                logging.warning("Using found fake array %s" % my_array[0])
            my_array = my_array[0]
    else:
        my_array = my_array[0]

    if my_array not in all_arrays:
        my_array = None
    elif all_arrays.index(my_array):
        all_arrays.remove(my_array)

    return my_array, all_arrays


def get_bam_sample_name(alignment_file):
    """Pull sample name from BAM header.

    Args:
        alignment_file (pysam.AlignmentFile):

    Returns:
        str
    """
    try:
        sample_name = alignment_file.header["RG"][0]["SM"]
        return sample_name
    except (KeyError, IndexError):
        try:  # just in case
            sample_name = re.search("\tSM:(\w*)\t?", alignment_file.text).groups()[0]
            return sample_name
        except IndexError:
            logging.error("Couldn't pull sample name from bam header")
            exit(1)


def main():
    args = parse_args()

    if args.user_input is not None:
        make_commands(args)
        exit(0)

    probe_sites = build_probe_sites(args.probelist)
    symlink_is_set = False
    mode = None

    if not args.in_silico:
        # Get the inputs
        if args.bam.endswith(".bam"):
            mode = "rb"
        elif args.bam.endswith(".cram"):
            mode = "rc"
        else:
            logging.error("Provided alignment file %s has unrecognized extension. Expecting .cram or .bam" % args.bam)
            exit(1)

        # There is a bug in pysam 0.8.4 for finding cram indexes, so we make a symlink and remove it at the end of run.
        # https://github.com/pysam-developers/pysam/pull/199/files#diff-b10e65edffc319543a90ba1ebf70e323

        if mode == "rc":
            index_path = None

            if os.path.exists(args.bam + '.crai'):
                index_path = args.bam + '.crai'
            elif os.path.exists(args.bam[:-1] + 'i'):
                index_path = args.bam[:-1] + 'i'
            elif os.path.exists(args.bam[:-4] + 'csi'):
                index_path = args.bam[:-4] + 'csi'
            else:
                logging.error("Provided alignment file %s doesn't have an index; expecting .cram or .bam" % args.bam)
                exit(1)

            if not os.path.exists(args.bam + '.crai'):
                symlink_is_set = True
                os.symlink(index_path, args.bam + '.crai')

            bam = pysam.AlignmentFile(args.bam, mode)

        else:
            bam = pysam.AlignmentFile(args.bam, mode)

        if not bam.has_index():
            logging.error("Provided alignment file %s doesn't have an index. Expecting .cram or .bam" % args.bam)
            exit(1)

        sample_name = get_bam_sample_name(bam)
        in_silico, bam_probe_report = probe_bam(bam, probe_sites)

        if args.genotype_output:
            save_in_silico_array(in_silico, bam_probe_report, sample_name, args.genotype_output)
            logging.info("Genotypes at probe sites in bam dumped to {}".format(args.genotype_output))
            exit(0)

    else:
        in_silico, bam_probe_report, sample_name = load_in_silico_array(args.in_silico)

    args.output.write("##ErisReport             %s\n" % sample_name)

    if args.Array.endswith('_FAKE.birdseed'):
        with open(args.Array, 'w') as array_file:
            logging.info("Writing fake array file '{}'".format(args.Array))
            array_file.write('')

    if args.Array == 'None':
        args.Array = None

    self_array_fn, all_array_fns = find_birdseed(sample_name, args.array_dir, args.Array)

    if self_array_fn is not None:
        logging.info("Found sample matching array %s" % self_array_fn)
        self_array = build_array(probe_sites, self_array_fn)
    else:
        logging.info("No sample matching array found (%s)" % sample_name)
        self_array = None

    if self_array is not None:
        contamin_report = generate_contamination_report(in_silico, self_array, probe_sites)
        self_concord = make_concordance(probe_sites, in_silico, build_array(probe_sites, self_array_fn))
    else:
        contamin_report = ContaminationReport()
        self_concord = ConcordanceReport()

    # output bam_probe_report
    args.output.write(bam_probe_report.pretty_report() + '\n')
    all_reports = {}
    if self_array is not None:
        # output contamination
        args.output.write(contamin_report.pretty_report() + '\n')
        # output self concordance report
        array_sample_name = os.path.basename(self_array_fn).replace('.birdseed', '')
        args.output.write("##SelfArray              %s\n" % array_sample_name)
        args.output.write(self_concord.pretty_report() + '\n')
        all_reports[array_sample_name] = self_concord
    else:
        args.output.write("##SelfArrayNotFound      null\n")
        # output no self found
        pass

    logging.info("Checking concordance against %d arrays" % (len(all_array_fns)))
    for fn in all_array_fns:
        array_sample_name = os.path.basename(fn).replace('.birdseed', '')
        array_sites = build_array(probe_sites, fn)
        all_reports[array_sample_name] = make_concordance(probe_sites, in_silico, array_sites)
        logging.debug("Finished checking " + fn)

    # sort all reports
    # Name, value
    sorted_reports = sorted(all_reports.items(), key=lambda y: y[1].concordance, reverse=True)
    judgement, avg_concordance, z_score = run_judgement(sorted_reports, self_concord, contamin_report)

    args.output.write("##AverageConcordance     %.2f%%\n" % avg_concordance)
    args.output.write("##SelfConcordanceZscore  %.2f\n" % z_score)
    args.output.write("##Judgement              %s\n" % judgement)

    # write our report
    # should sort all_reports based on Concordance.concordance and output
    args.output.write("#Sample\t" + sorted_reports[0][1].get_header() + '\n')
    for i in sorted_reports:
        args.output.write(i[0] + '\t' + str(i[1]) + '\n')

    args.output.close()

    if mode == 'rc' and os.path.exists(args.bam + '.crai') and symlink_is_set:
        logging.info("Found .cram.crai")
        if os.path.islink(args.bam + '.crai'):
            logging.info("Removing symlink")
            os.remove(args.bam + '.crai')

    # Removing fake arrays causes race condition with other eris jobs running
    # with the same array directory.
    # I don't see any harm in leaving the fake arrays where they are.

    # if args.Array.endswith('_FAKE.birdseed'):
        # logging.info("Removing fake array file '{}'".format(args.Array))
        # os.remove(args.Array)

    logging.info("Finished")


if __name__ == '__main__':
    main()
