#!/usr/bin/perl

########################
### Required Modules ###

use lib "/home/ya76/BioPerl/";
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

#######################################################################################################################################################################
#   
#   PAML_branchSite.pl 
#					
#	Requirements:
#			
#			1. In addition to the BioPerl modules required for this script to work (see code header), the following programs/scripts are required (make sure they are in your path variable):
#				
#				- ucsc utilities (faOneRecord & faSomeRecords, found here -> http://hgdownload.soe.ucsc.edu/admin/exe/)
#				- fastagrep.pl
#				- fastaSortByName.pl
#				- paml_prep.pl
#
#			2. CDS file should contain all sequences to be analyzed. Transcript ID within species (among lines) should be the same. Header should have follwoing format:
#				>transcript_id .......... species=<species_id>
#
#			
#   Required arguments:
#
#	--CDS_file|c <string>				sequence file containing CDS from multiple species/lines.
#
#	--gene_list|g <string>				list file of genes to be analyzed.
#		OR
#	--gene_id|i <string>				transcript/gene to be analyzed
#   
#	--output|o					name of directory to place outputs 
#
#
#
#   Options:
#
#	--show_samples					display species/line IDs in CDS_file and exit (only requires "-c <CDS_file>" argument).
#
#	--include_samples				include specified species/lines in --samples_file <string> (use line ID if present, otherwise use species ID).
#
#	--exclude_samples				exclude specified species/lines in --samples_file <string> (use line ID if present, otherwise use species ID).
#
#
#	--include_outgroups
#
#		--orthology_file <string>		file with orthology information:<transcript_id>\\t<orthologue_id>
#
#
#	--view_DNA_alignment				open DNA alignment in Geneious.
#
#	--view_Protein_alignment			open Protein alignment in Geneious.
#
#	--save_alignments				retain alignment files
#
#	--calculate_KaKs				Calculate pairwise Ka/Ks.
#
#	--run_PAML					Perform PAML codeml (will only perform branch-site test)
#
#		--ctl_dir				directory with control and tree files
#
#	--help						
#
#                           
#######################################################################################################################################################################


__EOUSAGE__

    ;


#####################################
## Set BioPerl alignment protocol ###

my $aln_factory = Bio::Tools::Run::Alignment::Clustalw->new;


###################################
### Define command line options ###

my $CDS_file;
my $gene_list;
my $gene_id;

my $output_dir;

my $show_samples_flag = 0;
my $include_samples_flag = 0;
my $exclude_samples_flag = 0;
my $samples_file;

my $include_outgroups_flag = 0;
my $orthology_file;

my $view_aln_dna_geneious_flag = 0;
my $view_aln_aa_geneious_flag = 0;
my $save_alignments_flag = 0;

my $calculate_KaKs_flag = 0;
my $run_PAML_flag = 0;
my $ctl_dir;

my $help_flag = 0;


&GetOptions (  

    'CDS_file|c=s' => \$CDS_file,
    'gene_list|g=s' => \$gene_list,
    'gene_id|i=s' => \$gene_id,
    'output|o=s' => \$output_dir,
    
    'show_samples' => \$show_samples_flag,
    'include_samples=s' => \$include_samples_flag,
    'exclude_samples=s' => \$exclude_samples_flag,
    'samples_file=s' => \$samples_file,
    
    'include_outgroups' => \$include_outgroups_flag,
    'orthology_file=s'=> \$orthology_file,
    'view_DNA_alignment' => \$view_aln_dna_geneious_flag,
    'view_Protein_alignment' => \$view_aln_aa_geneious_flag,
    'save_alignments' => \$save_alignments_flag,
    
    'calculate_KaKs' => \$calculate_KaKs_flag,
    'run_PAML' => \$run_PAML_flag,
    'ctl_dir=s' => \$ctl_dir,
    
    'help' => \$help_flag,
    
);


#####################################################
### Check command line arguments and housekeeping ###

# if command line argument are not recognized
if (@ARGV) {
    die "Error, don't understand parameters: @ARGV";
}

if ($help_flag) {
    die $usage;
}

# only show samples in CDS_file
if ($show_samples_flag) {
    system ("grep '>' $CDS_file | sed 's/.*species=/species=/g' | sort -u");
	exit(0);
}

# if Mandatory arguments not present
unless ($gene_list || $gene_id || $output_dir) {
	die "\nNo CDS file, gene ID or gene list file provided.\n\nAnd you need to specify output directory!\n $usage";
}

if (@ARGV) {
	die "Error, do not recognize params: @ARGV ";
}

# specify output dir cname
if ($output_dir) {
	mkdir($output_dir);
}

# if calculate KaKs, make tmp dir
if ($calculate_KaKs_flag) {
	system ("mkdir KaKs.tmp");
}

# if run PAML, make tmp dir
if ($run_PAML_flag) {
	unless ($ctl_dir) {
		die "\nError, Need to specify where PAML control file and tree file are, dude!";
	}
	system ("mkdir PAML.output");
}



######################################
###### Main execution of processes ###

main: {
				
		my @gene_object;
		
		###############################################################
		##### Set the gene object to a single gene or list of genes. ##
		if ($gene_id) {
			push (@gene_object, $gene_id);
		} elsif ($gene_list) {
			open (FILE1, $gene_list) or die ("Could not open file \n");
			while ($gene_list = <FILE1>){
				chomp $gene_list;
				my @file_bits = split(/\s+/,$gene_list);
				my $gene_names = $file_bits[0];
				push (@gene_object, $gene_names);
			}
		}
		
		
		###############################################################
		### Initiate the loop to analyze individual transcripts #######
	
		foreach my $transcript (@gene_object){
			
			###########################################################
			### Add orthologous sequences to the analysis...  #########
			### Processes FASTA header to include species and #########
			### line ID, if present.                          #########
			
			if ($include_outgroups_flag) {
				system ("grep $transcript $orthology_file | tr '\t' '\n' | sort -u > $transcript.list");
				system ("faSomeRecords $CDS_file $transcript.list $transcript.fa");
				
				#######################################################
				### Restrict to a subset of samples ###################
				
				if ($include_samples_flag) {
					system ("faOneRecord -f $samples_file $transcript.fa > $transcript.trimmed.fa");
         		}
         		
         		# header processing, contd.
				system ("sed -i.bak '/>/ s/ line=/_/g' $transcript.fa");
				system ("sed -i.bak '/>/ s/>.*species=/>/g' $transcript.fa");
				system ("rm $transcript.list *.bak");
### This part will be temporarily edited to alow MK test without adding 				
			##########################################################
			### Don't add orthologous seqnuences to analysis. ########
       	 	} else {
				system ("faOneRecord $CDS_file $transcript > $transcript.fa");
				system ("sed -i.bak '/>/ s/ line=/_/g' $transcript.fa");
				system ("sed -i.bak '/>/ s/>.*species=/>/' $transcript.fa");
			}
			
			
			##########################################################
			### BioPerl process for sequence analysis ################
			
			# read in sequence	
			my $seqio = Bio::SeqIO->new(-file => "$transcript.fa",
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
					warn("Oi! $transcript has a stop codon!!!");
				}
				
				# Tcoffee can't handle '*' even if it is trailing
				$pseq =~ s/\*//g;
				$protein->seq($pseq);
				push @prots, $protein;
			}
 			
 			# Warn if only 1 sequence present
			if( @prots < 2 ) {
				warn("$transcript - Need at least 2 CDS sequences to proceed");
			}

			# Align the sequences with clustalw
			my $aa_aln = $aln_factory->align(\@prots);
			
			# project the protein alignment back to CDS coordinates
			my $dna_aln = aa_to_dna_aln($aa_aln, \%seqs); 
			my @each = $dna_aln->each_seq();
			
			# output alignments for downstream analysis
			my $out_dna = Bio::AlignIO->new(-file => ">$transcript.aln_dna.afa" ,
										-format => 'fasta');
			$out_dna -> write_aln($dna_aln);

			my $out_aa = Bio::AlignIO->new(-file => ">$transcript.aln_aa.afa" ,
										-format => 'fasta');
			$out_aa -> write_aln($aa_aln);
			
			# clean-up alignment header files
			system ("sed -i.bak '/>/ s/\\/.*//g' $transcript.*.afa");
			system ("rm *bak");
			
			
			############################################################
			### Perform the analyses specified in the command line #####
			
			# Calculate Ka/Ks (Adapted from BioPerl)
			if ($calculate_KaKs_flag) {
				open(OUT, ">$transcript.KaKs.txt") ||  die("cannot open output align_output for writing");

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
						print OUT join("\t", $transcript,$otus[$i]->display_id,
							$otus[$j]->display_id,$MLmatrix->[$i]->[$j]->{'dN'},
							$MLmatrix->[$i]->[$j]->{'dS'},
							$MLmatrix->[$i]->[$j]->{'omega'},
							sprintf("%.2f",$sub_aa_aln->percentage_identity),
							sprintf("%.2f",$sub_dna_aln->percentage_identity),
							), "\n";
					}
				}
				system ("mv $transcript.KaKs.txt KaKs.tmp");
			}
###############################################################################################			
			# Run PAML analyses
			if ($run_PAML_flag) {
			system ("cp $ctl_dir/* .");
			system ("fastaSortByName.pl $transcript.aln_dna.afa > $transcript.dnaP.afa");
			system ("fastaSortByName.pl $transcript.aln_aa.afa > $transcript.aaP.afa");
			system ("sed -i.bak '/>/ s/\\/.*//g' $transcript.*P.afa");
			system ("paml_prep.pl $transcript.aaP.afa $transcript.dnaP.afa -nogap -output paml > $transcript.phy");
			system ("sed -i.bak 's/paml.phy/$transcript.phy/g' codeml.ctl");
			
			### 1. NULL MODEL
			system ("sed -i.bak 's/paml.out/$transcript.H0.out/g' codeml.ctl");
			system ("codeml");
			my $H0_omega = `grep omega $transcript.H0.out | awk '{print \$4}'`;
			my $H0_lnL = `grep -w -m 1 'lnL' $transcript.H0.out | awk '{print \$5}'`;
			system ("sed -i.bak 's/$transcript.H0.out/paml.out/g' codeml.ctl");
			system ("mv $transcript.H0.out PAML.output");
			
			### BRANCH-SITE MODELS
			system ("sed -i.bak 's/NSsites = 0/NSsites = 2/g' codeml.ctl");
			system ("sed -i.bak 's/model = 0/model = 2/g' codeml.ctl");
			system ("sed -i.bak 's/fix_omega = 0/fix_omega = 1/g' codeml.ctl");
				# 4a. Damr
				# edit tree
				system ("sed -i.bak 's/Damr /Damr \#1/g' ANLV.tree");
				# edit codeml.ctl 
				system ("sed -i.bak 's/paml.out/$transcript.brSt_Damr.H0.out/g' codeml.ctl");
				### run and get H0 info
				system ("codeml");
				my $Damr_brSt_H0 = `grep -w -m 1 'lnL' $transcript.brSt_Damr.H0.out | awk '{print \$5}'`;
				### swith to H1
				system ("sed -i.bak 's/$transcript.brSt_Damr.H0.out/$transcript.brSt_Damr.H1.out/g' codeml.ctl");
				system ("sed -i.bak 's/fix_omega = 1/fix_omega = 0/g' codeml.ctl");
				### run and get H1 info
				system ("codeml");
				my $Damr_brSt_H1 = `grep -w -m 1 'lnL' $transcript.brSt_Damr.H1.out | awk '{print \$5}'`;
				### back to normal
				system ("sed -i.bak 's/$transcript.brSt_Damr.H1.out/paml.out/g' codeml.ctl");
				system ("sed -i.bak 's/fix_omega = 0/fix_omega = 1/g' codeml.ctl");
				system ("sed -i.bak 's/Damr \#1/Damr /g' ANLV.tree");

				# 4b. Dnov
				# edit tree
				system ("sed -i.bak 's/Dnov /Dnov \#1/g' ANLV.tree");
				# edit codeml.ctl 
				system ("sed -i.bak 's/paml.out/$transcript.brSt_Dnov.H0.out/g' codeml.ctl");
				### run and get H0 info
				system ("codeml");
				my $Dnov_brSt_H0 = `grep -w -m 1 'lnL' $transcript.brSt_Dnov.H0.out | awk '{print \$5}'`;
				### swith to H1
				system ("sed -i.bak 's/$transcript.brSt_Dnov.H0.out/$transcript.brSt_Dnov.H1.out/g' codeml.ctl");
				system ("sed -i.bak 's/fix_omega = 1/fix_omega = 0/g' codeml.ctl");
				### run and get H1 info
				system ("codeml");
				my $Dnov_brSt_H1 = `grep -w -m 1 'lnL' $transcript.brSt_Dnov.H1.out | awk '{print \$5}'`;
				### back top normal
				system ("sed -i.bak 's/$transcript.brSt_Dnov.H1.out/paml.out/g' codeml.ctl");
				system ("sed -i.bak 's/fix_omega = 0/fix_omega = 1/g' codeml.ctl");
				system ("sed -i.bak 's/Dnov \#1/Dnov /g' ANLV.tree");
				
				# 4c. Dvir
				# edit tree
				system ("sed -i.bak 's/Dvir /Dvir \#1/g' ANLV.tree");
				# edit codeml.ctl 
				system ("sed -i.bak 's/paml.out/$transcript.brSt_Dvir.H0.out/g' codeml.ctl");
				### run and get H0 info
				system ("codeml");
				my $Dvir_brSt_H0 = `grep -w -m 1 'lnL' $transcript.brSt_Dvir.H0.out | awk '{print \$5}'`;
				### swith to H1
				system ("sed -i.bak 's/$transcript.brSt_Dvir.H0.out/$transcript.brSt_Dvir.H1.out/g' codeml.ctl");
				system ("sed -i.bak 's/fix_omega = 1/fix_omega = 0/g' codeml.ctl");
				### run and get H1 info
				system ("codeml");
				my $Dvir_brSt_H1 = `grep -w -m 1 'lnL' $transcript.brSt_Dvir.H1.out | awk '{print \$5}'`;
				### back top normal
				system ("sed -i.bak 's/$transcript.brSt_Dvir.H1.out/paml.out/g' codeml.ctl");
				system ("sed -i.bak 's/Dvir \#1/Dvir /g' ANLV.tree");
				
				# Damr-Dnov
				# edit tree
				system ("sed -i.bak 's/Dnov )/Dnov )\$1/g' ANLV.tree");
				# edit codeml.ctl 
				system ("sed -i.bak 's/paml.out/$transcript.brSt_Damr-Dnov.H0.out/g' codeml.ctl");
				### run and get H0 info
				system ("codeml");
				my $DamrNov_brSt_H0 = `grep -w -m 1 'lnL' $transcript.brSt_Damr-Dnov.H0.out | awk '{print \$5}'`;
				### swith to H1
				system ("sed -i.bak 's/$transcript.brSt_Damr-Dnov.H0.out/$transcript.brSt_Damr-Dnov.H1.out/g' codeml.ctl");
				system ("sed -i.bak 's/fix_omega = 1/fix_omega = 0/g' codeml.ctl");
				### run and get H1 info
				system ("codeml");
				my $DamrNov_brSt_H1 = `grep -w -m 1 'lnL' $transcript.brSt_Damr-Dnov.H1.out | awk '{print \$5}'`;
				### back top normal
				system ("sed -i.bak 's/$transcript.brSt_Damr-Dnov.H1.out/paml.out/g' codeml.ctl");
				system ("sed -i.bak 's/Dnov )\$1/Dnov )/g' ANLV.tree");
		
			system ("sed -i.bak 's/NSsites = 2/NSsites = 0/g' codeml.ctl");
			system ("sed -i.bak 's/model = 2/model = 0/g' codeml.ctl");
			system ("mv $transcript.brSt_*.out PAML.output");
			system ("sed -i.bak 's/$transcript.phy/paml.phy/g' codeml.ctl");	
			
			
			################## COMBINE OUTPUTS ######################
			my $outfile = "$transcript.branchSite.tmp.txt";
			my $out;
			open($out, ">$outfile") or die "can't open $outfile: $!\n";
			print $out "$transcript\t$H0_omega\t$H0_lnL\t$Damr_brSt_H0\t$Damr_brSt_H1\t$Dnov_brSt_H0\t$Dnov_brSt_H1\t$Dvir_brSt_H0\t$Dvir_brSt_H1\t$DamrNov_brSt_H0\t$DamrNov_brSt_H1\n";
		

 			system ("rm $transcript.*P.afa");
 			system ("rm $transcript.phy");
			system ("cat $transcript.branchSite.tmp.txt | awk -v OFS='\t' '\$1=\$1' | tr '\n' '\t' | perl -p -e 's/\t\$/\n/g' > PAML.output/$transcript.branchSite.txt");
			system ("rm $transcript.branchSite.tmp.txt");
			system ("rm ANLV.tree codeml.ctl");
			}
			
				
			
			
			# option to save the alignment files
			if ($save_alignments_flag || $view_aln_dna_geneious_flag || $view_aln_aa_geneious_flag) {
				system ("mv $transcript.aln_dna.afa $output_dir/$transcript.DNA_alignment.afa");
				system ("mv $transcript.aln_aa.afa $output_dir/$transcript.Protein_alignment.afa");
			} else {
				system ("rm $transcript.aln*.afa");
			}
			
			# option to view alignment files in Geneious
			if ($view_aln_dna_geneious_flag) {
				system ("open -a Geneious $output_dir/$transcript.DNA_alignment.afa");
			}
				
			
			if ($view_aln_aa_geneious_flag) {
				system ("open -a Geneious $output_dir/$transcript.Protein_alignment.afa");
			}
			
			# clean-up
			system ("rm $transcript.fa");
			
			##########################################################
			############### END OF LOOP ##############################
		}		
		
		# Format final output files
		
		if ($calculate_KaKs_flag) {
			system ("cat KaKs.tmp/*KaKs.txt | sed '/TRANSCRIP/d' > KaKs.tmp/All.tmp.KaKs.txt");
			system ("echo 'TRANSCRIPT\tSEQ1\tSEQ2\tKa\tKs\tKa/Ks\tPROT_PERCENTID\tCDNA_PERCENTID' | cat - KaKs.tmp/All.tmp.KaKs.txt > $output_dir/FINAL.KaKs.txt");
			system ("rm -r KaKs.tmp");
		}
		
		if ($run_PAML_flag){
				system ("echo 'TRANSCRIPT\tomega\tH0_lnL\tDamr_brSt_H0\tDamr_brSt_H1\tDnov_brSt_H0\tDnov_brSt_H1\tDvir_brSt_H0\tDvir_brSt_H1\tDamrNov_brSt_H0\tDamrNov_brSt_H1' | cat - PAML.output/*.branchSite.txt > PAML.output/branchSite.all.txt");
				system ("sed -i.bak '/^\$/d' PAML.output/branchSite.all.txt");
				system ("rm *bak PAML.output/*bak PAML.output/*.branchSite.txt");
				system ("mv PAML.output $output_dir/");
				system ("rm r* 2N* 4* lnf");
		}	
}
