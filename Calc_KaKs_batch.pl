#!/usr/bin/perl

=head1 SYNOPSIS

usage: Calc_KaKs_batch.pl <transcript_list> <alignment_file> <sample_pair> <header>


 	- transcript_list: list of transcripts to calculate
 	- alignment_file: file containing all alignments in FASTA format
 	- sample pair: file containing three-letter code for species name, each in newline.
 	
KaKs_calculator parameter is set to "Model Averaging on a set of candidate models" (MA). Edit the original script file to choose different method.


  
=cut

use strict;
use warnings; 
use Pod::Usage;

my $input_file = $ARGV[0];
my $fasta_alignment_file = $ARGV[1];
my $sample_pair = $ARGV[2];
my $header = $ARGV[3];
my ($help);


pod2usage(1) if($help);
pod2usage("No files given!\n")  if ((!$input_file || !$fasta_alignment_file));



open (FILE1, $input_file) or die ("Could not open file \n");

while ($input_file = <FILE1>){
chomp $input_file;
my @file_bits = split(/\s+/,$input_file);
my $dir_name = $file_bits[0];

my @array = ();
push (@array, $dir_name);

foreach my $int (@array){

########## USE this script to get the alignment for PAML ######
######### Sequences must be pre-aligned #######

system ("faOneRecord $fasta_alignment_file $int > $int.tmp.fa");
# system ("fastaSortByName.pl $int.tmp.fa > $int.fa");

system ("sed -i.bak '/amr/ s/>/>amr/' $int.tmp.fa");
system ("sed -i.bak '/lum/ s/>/>lum/' $int.tmp.fa");
system ("sed -i.bak '/nov/ s/>/>nov/' $int.tmp.fa");
system ("sed -i.bak '/vir/ s/>/>vir/' $int.tmp.fa");

system ("sed -i.bak 's/TCONS.*//g' $int.tmp.fa");

system ("faSomeRecords $int.tmp.fa $sample_pair $int.fa");

system ("parseFastaIntoAXT.pl $int.fa");
system ("KaKs_Calculator -i $int.fa.axt -o $int.fa.axt.kaks");
my $KaKs_row = `grep -v 'Sequence' $int.fa.axt.kaks`;

my $outfile = "$int.KaKs.txt";
my $out;
open($out, ">$outfile") or die "can't open $outfile: $!\n";

print $out "$int\t$KaKs_row";

# system ("rm $int.tmp.fa.bak");
# system ("rm $int.fa");
# system ("rm $int.tmp.fa");
# system ("rm $int.fa.axt");
# system ("rm $int.fa.axt.kaks");

			}


	}
close FILE1;

