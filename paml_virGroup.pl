#!/usr/bin/perl -w

=head1 SYNOPSIS

 paml_virGroup.pl <gene_list> <alignment_file> <codeml_ctrl_file> <out_prefix> <branch_in> <branch_out> <header_file>

  
=cut

use strict;
use warnings; 
use Bio::Align::Utilities qw(aa_to_dna_aln);
use Bio::SeqIO;
use Bio::AlignIO;
use Pod::Usage;

my $input_file = $ARGV[0];
my $fasta_alignment_file = $ARGV[1];
my $codeml_ctrl_file = $ARGV[2];
my $branch_in = $ARGV[4];
my $branch_out = $ARGV[5];
my $header = $ARGV[6];
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
system ("fastaSortByName.pl $int.tmp.fa > $int.fa");

my $seqin = Bio::SeqIO->new(-file => "$int.fa",
                            -format => "fasta");

my $seqout = Bio::SeqIO->new(-file => ">$int.FASTA",
                            -format => "fasta");

# process each sequence
while ( my $seq = $seqin->next_seq ) {
    # translate them into protein
    my $protein = $seq->translate();
    $seqout->write_seq($protein);

	}
system ("sed -i '/amr/ s/>/>amr/' $int.fa");
system ("sed -i '/lum/ s/>/>lum/' $int.fa");
system ("sed -i '/nov/ s/>/>nov/' $int.fa");
system ("sed -i '/vir/ s/>/>vir/' $int.fa");

system ("sed -i '/amr/ s/>/>amr/' $int.FASTA");
system ("sed -i '/lum/ s/>/>lum/' $int.FASTA");
system ("sed -i '/nov/ s/>/>nov/' $int.FASTA");
system ("sed -i '/vir/ s/>/>vir/' $int.FASTA");

system ("sed -i 's/TCONS/ TCONS/g' $int.fa");
system ("sed -i 's/TCONS/ TCONS/g' $int.FASTA");

system ("paml_prep.pl $int.FASTA $int.fa -nogap -output paml > $int.phy");
system ("perl -pi -e 's/paml.phy/$int.phy/g' $codeml_ctrl_file");
system ("perl -pi -e 's/paml.out/$int.out/g' $codeml_ctrl_file");
system ("codeml $codeml_ctrl_file");
system ("perl -pi -e 's/$int.phy/paml.phy/g' $codeml_ctrl_file");
system ("perl -pi -e 's/$int.out/paml.out/g' $codeml_ctrl_file");

#### Modify this line to grep your desired output from the PAML *.out file ####

my $log_lik = `grep -w -m 1 'lnL' $int.out | awk '{print \$5}'`;
my $val = `grep -w '$branch_in' $int.out | grep -vw '$branch_out' | tr '\n' ' '`;
 

###############################################################################

my $outfile = "$int.dNdS.txt";
my $out;
open($out, ">$outfile") or die "can't open $outfile: $!\n";

print $out "$int\t$val\t$log_lik\n";
system ("rm $int.tmp.fa");
system ("rm $int.fa");
system ("rm $int.FASTA");
system ("rm $int.out");
system ("rm $int.phy");
			}


	}
my $prefix = $ARGV[3];
system ("cat *dNdS.txt > $prefix.dNdS.all.tmp1.txt");
system ("awk -v OFS='\t' '\$1=\$1' $prefix.dNdS.all.tmp1.txt > $prefix.dNdS.all.tmp.txt");
system ("rm *dNdS.txt *tmp1.txt");
system ("cat $header $prefix.dNdS.all.tmp.txt > $prefix.dNdS.all.txt");
system ("sed -i '/^\$/d' $prefix.dNdS.all.txt");
system ("rm $prefix.dNdS.all.tmp.txt");
system ("rm r* 2N* 4* lnf");

close FILE1;

