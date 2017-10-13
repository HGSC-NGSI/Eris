#! /hgsc_software/perl/latest/bin/perl


#############################################
#
#	Create cSNP birdseed from .csv file
#
#
#
#	Created by:David Sexton
#
#############################################

use strict;
use warnings;

open (FIN_SNP, "< $ARGV[0]") or die "No array file path entered: $!\n";

#open (FOUT_SNP, "> $ARGV[1]") or die "No birdseed file path entered: $!\n"; 

open (FIN_PROBE, "< $ARGV[1]") or die "No probelist path entered: $!\n";

my %pr_data;

while (<FIN_PROBE>) {

    my $line = $_;
	
    chomp $line;
	    
    my @fields = split(/\t/, $line);
		
    my ($rs_id, $chromosome, $position, $ref_allele, $var_allele) = ($fields[2], $fields[0], $fields[1], $fields[5], $fields[6]);
		    
    push (my @alleles, ($chromosome, $position, $ref_allele, $var_allele));
			
    push( @{$pr_data{$rs_id}}, @alleles);

}

while (<FIN_SNP>) {
	
	chomp $_;
	
	my @fields = split (/[,\t]/, $_); #changed to read csv and txt files
	
	my $field_size = scalar (@fields);
	
	if ($field_size >= 10) {
	
		my ($sample_id, $rs_id, $chromosome, $position, $allele1, $allele2) = ($fields[0],$fields[1],$fields[2], $fields[3], $fields[5], $fields[6]);
		
		$sample_id = $sample_id.".birdseed";
		
		$rs_id =~ s/(exm-|exm)//g;
		
		if ($rs_id =~ m/IND/g || $sample_id =~ m/Sample/g || $allele1 =~ m/-/g || $allele2 =~ m/-/g) {
			
			next;
			
		}elsif (defined @{$pr_data{$rs_id}}) {
			
		    #print "working $rs_id\n";
	
		    #this means that data is appending to the file, so you will need to delete file between runs.

		    open (FOUT_SNP, ">> $sample_id") or die "cannot open birdseed output: $!\n";
		    	
		    my @alleles = @{$pr_data{$rs_id}};
			
			#$refvar =~ m/\[(\w)\/(\w)\]/g;
	
			#my $ref = $1;
			
			#my $var = $2;
			
			#    print "$allele1\t$allele2\t$alleles[3]\n";
			if ($allele1 =~ m/$alleles[2]/g || $allele2 =~ m/$alleles[2]/g || $allele1 =~ m/$alleles[3]/g || $allele2 =~ m/$alleles[3]/g) {
			   
			    my $genotype = $allele1.$allele2;
	
			    print FOUT_SNP "$chromosome\t$position\t$alleles[2]\t$genotype\n";

			}else {
				
				$allele1 =~ tr/ACTG/TGAC/;
				
				$allele2 =~ tr/ACTG/TGAC/;
				
				my $genotype = $allele1.$allele2;
				
				print FOUT_SNP "$chromosome\t$position\t$alleles[2]\t$genotype\n";
				
			}
		    close FOUT_SNP;
	
		}	

	}else {
		
		next;
		
	}
	
}

close FIN_SNP;

#close FOUT_SNP;

close FIN_PROBE;
