# Eris

Eris is a program for calculating sample concordance between BAM/CRAM
alignment files and SNP array files in birdseed format. Also provided are
scripts for converting array files for some common SNP chips (such as
Fluidigm SNPTrace and OmniExpress1.0) into the birdseed format.

The output report files provide an estimated level of contamination
between a given BAM and its specific array file and a sorted table of
estimated concordance for each array file in a given directory.

Eris can also be run with the -g flag to create a file containing
genotypes at each probelist site from the given BAM. These files are
much smaller than the BAM itself and can be used with the -i flag to run
concordance in place of the BAM.

Contributors:
    Adam English, Adam Mansfield, Shruthi Ambreth, Jeffrey Reid

### Requirements

Requires Python 3 and pysam 8.4+.

Some of the conversion scripts require Perl 5.

### Usage
* Use provided scripts (under array_kits/) to convert array files
into birdseed format.

* Move all birdseed files into one directory.

* Ensure the BAM/CRAM file you will use is indexed (such as with
samtools index).

* Then run:
```
python3 Eris.py -b BAM -a ARRAY_DIR -A SELF_ARRAY -p PROBELIST -o OUTFILE
```

* For example:
```
python3 Eris.py -b NA12878.bam -a /path/to/birdseeds_dir \
      -A /path/to/birdseeds_dir/NA12878.birdseed -o NA12878.report.txt \
      -p ../array_kits/fluidigm_37/probelist.txt
```

### Copyright

Copyright 2017 Baylor College of Medicine Human Genome Sequencing Center
