#!/usr/bin/env perl

#Usage: make_phylip_trees.pl <gene_list.file> <fasta_all_alignments> <PHYLIP_parameter_file>


use strict;
use warnings;

my $input_file = $ARGV[0];
my $fasta_alignment_file = $ARGV[1];
my $param_file = $ARGV[2];

open (FILE1, $input_file) or die ("Could not open file \n");

while ($input_file = <FILE1>){
chomp $input_file;
my @file_bits = split(/\s+/,$input_file);
my $dir_name = $file_bits[0];

my @array = ();
push (@array, $dir_name);

foreach my $int (@array){

########## ######
######### fetch sequences, edit accession, convert to phylip, then make tree #######

system ("faOneRecord $fasta_alignment_file $int > $int.fa");

system ("sed -i '/amr/ s/>/>amr/' $int.fa");
system ("sed -i '/lum/ s/>/>lum/' $int.fa");
system ("sed -i '/nov/ s/>/>nov/' $int.fa");
system ("sed -i '/vir/ s/>/>vir/' $int.fa");

system ("fasta_to_phylip.sh $int.fa > $int.phy");
system ("mv $int.phy infile");

system ("dnaml < $param_file > $int.terminal_file");
system ("mv outtree $int.outtree");
system ("mv outfile $int.outfile");


system ("rm $int.fa");
system ("rm infile");
			}

	}

close FILE1;

my $prefix = $ARGV[0];
system ("cat *.outtree > $prefix.ALL.trees.newick");
system ("mkdir $prefix.Make.Trees.OutFiles");
system ("mkdir $prefix.Make.Trees.OutFiles/tree.files.by.gene");
system ("mkdir $prefix.Make.Trees.OutFiles/out.files.by.gene");
system ("mkdir $prefix.Make.Trees.OutFiles/terminal.files.by.gene");
system ("mv *.outfile $prefix.Make.Trees.OutFiles/out.files.by.gene");
system ("mv *.outtree $prefix.Make.Trees.OutFiles/tree.files.by.gene");
system ("mv *.terminal_file $prefix.Make.Trees.OutFiles/terminal.files.by.gene");


