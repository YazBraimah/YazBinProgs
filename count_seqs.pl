#!/usr/bin/perl

=head1 SYNOPSIS

count how many sequences per transcript:

count_seqs.pl <gene_list> <cds_file>

  
=cut

use strict;
use warnings;
use Pod::Usage;

my $gene_list = $ARGV[0];
my $cds_file = $ARGV[1];
my ($help);

pod2usage(1) if($help);
pod2usage("No files given!\n")  if ((!$gene_list || !$cds_file));

open (FILE1, $gene_list) or die ("Could not open file \n");

system ("mkdir tmp.numbers");

while ($gene_list = <FILE1>){
chomp $gene_list;
my @file_bits = split(/\s+/,$gene_list);
my $dir_name = $file_bits[0];

my @array = ();
push (@array, $dir_name);

foreach my $int (@array){

my $seq_num = `faOneRecord $cds_file $int | grep '>' | wc -l`;

my $outfile = "$int.number.txt";
my $out;
open($out, ">$outfile") or die "can't open $outfile: $!\n";

print $out "$int\t$seq_num";

system ("mv $int.number.txt tmp.numbers");

	}
	
}

system ("cat tmp.numbers/*.number.txt > all.transcript.numbers.txt");
system ("rm -r tmp.numbers");

close FILE1;