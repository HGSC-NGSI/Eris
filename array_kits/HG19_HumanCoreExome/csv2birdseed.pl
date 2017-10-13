#!/usr/bin/perl
use strict;
#use warnings;

if(scalar @ARGV != 2){
	print "perl csnpcsv2birdseed.pl /path/to/probelist /path/to/csv/file\n";
	exit(0);
}
my $probelist = $ARGV[0];
my $raw_bs = $ARGV[1]; 

sub rev_complement{
        my $string = shift;
        $string =~ tr/ATGC/TACG/;
        my $rev_compl = scalar reverse("$string");
        return $rev_compl;
}

print STDOUT "Reading $probelist\n";
my %probes;
open(FIN,"<$probelist") or die $!;
while(<FIN>){
   #1       60520956        64173   GGTTGTTCTCTGGGT CTTTTTCTCTTTCAC G       A       0.998   0.001   0
   chomp;
   my @fields = split ("\t",$_);
   my $chr = $fields[0];
   $chr =~ s/[Cc]hr//;
   my $pos = $fields[1];
   $probes{$chr."_".$pos} = $_;
}
close(FIN);
print STDOUT "Probelist size ". scalar keys %probes;
print STDOUT "\n";


#Sample ID,SNP Name,Chr,Position,Allele1 - Design,Allele1 - Top,Allele2 - Top,SNP,Allele1 - Plus,Allele2 - Plus,GC Score,B Allele Freq,Log R Ratio
#6142-PT01-N-1_2,exm1090324,14,23848688,C,G,G,[T/C],C,C,0.5273,0.9907,-0.0848
print STDOUT "Reading raw birdseed $raw_bs\n";
my $birdseed;
my $count  = 0;
open(FINR,"<$raw_bs") or die $!;
while(<FINR>){
   chomp;
   my @fields = split(/[,\t]/,$_);
   if($fields[7] =~ m/\[\w+\/\w+\]/){
	
	$fields[0] !~ s/([^\d+\w+\-\_\s+])//g; #remove characters that are not white space/digits/character/underscore/hyphen
	my @birdseed_name = split(" ",$fields[0]);
	$birdseed = join("_",@birdseed_name).".birdseed";
	open(FOUT,">>$birdseed");

	my $rsid = $fields[1];
	my $chr = $fields[2];
	my $pos = $fields[3];
	my $allele1 = $fields[8];
	my $allele2 = $fields[9];

	$rsid =~ s/(exm-|exm)//g;
	$chr =~ s/[Cc]hr//;
 	if(exists $probes{$chr."_".$pos}){
		my  ($chr_v,$pos_v,$rsid_v,$ls,$rs,$ref,$var) = split("\t",$probes{$chr."_".$pos});
		if($allele1 =~ m/[$ref$var]/ and $allele2 =~ m/[$ref$var]/){	
			print FOUT "$chr\t$pos\t$ref\t$allele1$allele2\n";
		}
		else{
			print STDOUT "$chr:$pos $allele1/$allele2 $ref/$var\n";
		}
	}	
	close(FOUT);
   }
}
close(FINR);
