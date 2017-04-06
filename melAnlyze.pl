#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config no_ignore_case bundling pass_through);
use FindBin;
use File::Basename;
use Bio::Tools::Run::Phylo::PAML::Codeml;
use Bio::Tools::Run::Alignment::Clustalw;
# for projecting alignments from protein to R/DNA space
use Bio::Align::Utilities qw(aa_to_dna_aln);
# for input of the sequence data
use Bio::SeqIO;
use Bio::AlignIO;

my $usage = <<__EOUSAGE__;

#################################################################################### 
#   
#   seqAnlyze.pl 
#   
#  --CDS_file <string>					Sequence file containing CDS from multiple lines.
#
#  --gene_list <string>					List of genes to be analyzed.
#   
#
#
#	Options:
#
#  --view_DNA_alignment					Open DNA alignment in Geneious.
#
#  --view_Protein_alignment				Open DNA alignment in Geneious.
#
#  --calculate_KaKs						Calculate pairwise Ka/Ks.
#                           
####################################################################################



__EOUSAGE__

    ;

####### SET ALIGNMENT PROTOCOL ################
my $aln_factory = Bio::Tools::Run::Alignment::Clustalw->new;

my $CDS_file;
# my $annotation_file;
my $gene_list;

my $view_aln_dna_geneious_flag = 0;
my $view_aln_aa_geneious_flag = 0;

my $calculate_KaKs_flag = 0;

my $help_flag = 0;

&GetOptions (  
    
    ## Inputs and outputs
    'CDS_file|c=s' => \$CDS_file,
#     'annotation_file|a=s' => \$annotation_file, 
#     'output|o=s' => \$output_prefix,
    'gene_list|g=s' => \$gene_list,
    
    'view_DNA_alignment' => \$view_aln_dna_geneious_flag,
    'view_Protein_alignment' => \$view_aln_aa_geneious_flag,
    
    'calculate_KaKs' => \$calculate_KaKs_flag,
);

if (@ARGV) {
    die "Error, don't understand parameters: @ARGV";
}

if ($help_flag) {
    die $usage;
}

unless ($gene_list) {
    die $usage;
}

if (@ARGV) {
    die "Error, do not recognize params: @ARGV ";
}

main: {
		####### SETUP THE FOR-LOOP OF GENE LIST #######
		open (FILE1, $gene_list) or die ("Could not open file \n");
		while ($gene_list = <FILE1>){
			chomp $gene_list;
			my @file_bits = split(/\s+/,$gene_list);
			my $dir_name = $file_bits[0];
			my @array = ();
			push (@array, $dir_name);
			foreach my $int (@array){
				
				######### FETCH SEQUENCES AND EDIT SEQ HEADERS ###############
				system ("faOneRecord $CDS_file $int > $int.fa");
				# change sequence identifier to line identifier
				system ("sed -i.bak '/>/ s/FBtr.*line=//' $int.fa");
# 				system ("rm $int.fa.bak");
				
				######### CHECK THE SEQUENCES AND
				my $seqio = Bio::SeqIO->new(-file => "$int.fa",
										   -format => 'fasta');
				my %seqs;
				my @prots;
				# process each sequence
				while ( my $seq = $seqio->next_seq ) {
						$seqs{$seq->display_id} = $seq;
						# translate them into protein
						my $protein = $seq->translate();
						my $pseq = $protein->seq();
						if( $pseq =~ /\*/ &&
						$pseq !~ /\*$/ ) {
							 warn("$int - provided a CDS sequence with a stop codon.... what up with that?");
								# exit(0);
							}
						# Tcoffee can't handle '*' even if it is trailing
						$pseq =~ s/\*//g;
						$protein->seq($pseq);
						push @prots, $protein;
					}
 
				if( @prots < 2 ) {
						warn("$int - Need at least 2 CDS sequences to proceed");
					 # exit(0);
					}

				# Align the sequences with clustalw
				my $aa_aln = $aln_factory->align(\@prots);
				# project the protein alignment back to CDS coordinates
				my $dna_aln = aa_to_dna_aln($aa_aln, \%seqs); 
				my @each = $dna_aln->each_seq();

				########################################################################################################
				###################### THIS SEGMENT TO WRITE OUT SEQS ##################################################

				my $out_dna = Bio::AlignIO->new(-file => ">$int.aln_dna.tmp.afa" ,
											-format => 'fasta');
				$out_dna -> write_aln($dna_aln);

				my $out_aa = Bio::AlignIO->new(-file => ">$int.aln_aa.tmp.afa" ,
											-format => 'fasta');
				$out_aa -> write_aln($aa_aln);
				
				system ("sed -i.bak '/>/ s/\\/.*//g' $int.*.afa");
				
				if ($view_aln_dna_geneious_flag) {
						system ("open -a Geneious $int.aln_dna.tmp.afa");
						}
						
				if ($view_aln_aa_geneious_flag) {
						system ("open -a Geneious $int.aln_aa.tmp.afa");
						}				
				
				if ($calculate_KaKs_flag) {
						open(OUT, ">$int.KaKs.txt") ||  die("cannot open output align_output for writing");

						my $kaks_factory = Bio::Tools::Run::Phylo::PAML::Codeml->new
										( -params => { 'runmode' => -2,
														 'seqtype' => 1,
														 } );
 
						# set the alignment object
						$kaks_factory->alignment($dna_aln);
 
						# run the KaKs analysis
						my ($rc,$parser) = $kaks_factory->run();
						my $result = $parser->next_result;
						my $MLmatrix = $result->get_MLmatrix();
 
						my @otus = $result->get_seqs();
						# this gives us a mapping from the PAML order of sequences back to
						# the input order (since names get truncated)
						my @pos = map {
								my $c= 1;
								foreach my $s ( @each ) {
								last if( $s->display_id eq $_->display_id );
											$c++;
											}
											$c;
								} @otus;
 
						print OUT join("\t", qw(TRANSCRIPT SEQ1 SEQ2 Ka Ks Ka/Ks PROT_PERCENTID CDNA_PERCENTID)),"\n";
						foreach my $i ( 0 .. $#otus -1 ) {
							foreach my $j ( $i+1 .. $#otus ) {
							my $sub_aa_aln  = $aa_aln->select_noncont($pos[$i],$pos[$j]);
								my $sub_dna_aln = $dna_aln->select_noncont($pos[$i],$pos[$j]);
								print OUT join("\t", $int,$otus[$i]->display_id,
												$otus[$j]->display_id,$MLmatrix->[$i]->[$j]->{'dN'},
												$MLmatrix->[$i]->[$j]->{'dS'},
												$MLmatrix->[$i]->[$j]->{'omega'},
												sprintf("%.2f",$sub_aa_aln->percentage_identity),
												sprintf("%.2f",$sub_dna_aln->percentage_identity),
												), "\n";
												}
											}
						}
}
}
}