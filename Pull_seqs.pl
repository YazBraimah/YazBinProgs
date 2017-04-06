#!/software/perl/5.18.2/bin/perl

=head1 SYNOPSIS

Pull_seqs.pl <line_IDs_list_file> <gene_GTF_file> <transcripts_GTF_file>
  
=cut

use strict;
use warnings;
use Pod::Usage;
 
####### DEFINE COMMAND LINE ARGUMENTS #########
my $line_file = $ARGV[0];
my $genes_gtf_file = $ARGV[1];
my $transcripts_gtf_file = $ARGV[2];
my ($help);




####### SET HELP OUTPUT #######################
pod2usage(1) if($help);
pod2usage("No files given!\n")  if ((!$line_file || !$gtf_file));

####### create temporary directories
system ("mkdir CDS Proteins Transcripts Exons");

####### SETUP THE FOR-LOOP OF GENE LIST #######
open (FILE1, $line_file) or die ("Could not open file \n");
while ($line_file = <FILE1>){
chomp $line_file;
my @file_bits = split(/\s+/,$line_file);
my $dir_name = $file_bits[0];
my @array = ();
push (@array, $dir_name);
foreach my $int (@array){

	######### Generate CDS File ###############
	######### EDIT THIS SECTION FOR USE WITH OTHER SPECIES #######
	system ("gffread $genes_gtf_file -g $int.fa -w $int.exons.fa -x $int.CDS.fa -y $int.PEP.fa");
	system ("gffread $transcripts_gtf_file -g $int.fa -x $int.transcripts.fa");
	# Add strain identifier to seq files 
	system ("sed -i '/>/ s/\$/ line=$int/g' $int.*.fa");
	          

	system ("mv $int.exons.fa Exons");
	system ("mv $int.CDS.fa CDS");
	system ("mv $int.transcripts.fa Transcripts");
	system ("mv $int.PEP.fa Proteins");					
		}
	}



close FILE1;