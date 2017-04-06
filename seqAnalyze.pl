#!/usr/bin/env perl

########################
### Required Modules ###

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
#   seqAnlyze.pl 
#					
#	Population genetic analysis of coding sequences.
#
#		Requirements:
#			
#			1. In addition to the perl modules required for this script to work (see code header), the following programs/scripts are required:
#				
#				- ucsc utilities (faOneRecord & faSomeRecords, found here -> http://hgdownload.soe.ucsc.edu/admin/exe/)
#				- Perl scripts:
#
#					fastagrep.pl
#					MK.pl	
#					MK_trim.pl	
#
#			2. CDS file should contain all sequences to be analyzed. Transcript ID within species (among lines) should be the same. Header should have follwoing format:
#				>transcript_id .......... species=<species_id> line=<line_id>
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
#
#	--MK_test
#
#		--pol <string>				[pol|unpol] specify polarized or unpolarized analysis.
#		--outgroup <string>			[speciesID] 1 for unpolarized; must have 2 for polarized analysis.
#							(specificy two outgroups like so: --outgroup1 <string> --outgroup2 <string>)
#		--ingroup <string>			[speciesID] 1 species
#
#
#	--calculate_KaKs				Calculate pairwise Ka/Ks.
#
#	--run_PAML					Perform PAML codeml (for virilis group alignments only)
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

my $mk_test_flag = 0;
my $pol;
my $outgroup;
my $outgroup1;
my $outgroup2;
my $ingroup;

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
    
    'MK_test' => \$mk_test_flag,
    'pol=s' => \$pol,
    'outgroup=s' => \$outgroup,
    'outgroup1=s' => \$outgroup1,
    'outgroup2=s' => \$outgroup2,
    'ingroup=s' => \$ingroup,
    
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

# if MK test, check polarization options given and make tmp dir
if ($mk_test_flag) {
	unless ($pol =~ /^(pol|unpol)$/) {
		die "\nError, do not recognize polarization option [$pol]\nShould be pol or unpol";
	}
	system ("mkdir MK.tmp");	
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
					system ("fastagrep.pl -f $samples_file $transcript.fa > $transcript.trimmed.fa");
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
			
			# MK test
			if ($mk_test_flag) {
			
				# Edit/move alignment file to MK tmp dir
				system ("sed '/>/ s/_/; /g' $transcript.aln_dna.afa > MK.tmp/$transcript.fa ");
				system ("sed -i.bak '/>/ s/\$/;/g' MK.tmp/$transcript.fa ");
				
				# run MK.pl with "pol" option
				if ($pol eq "pol") {
					system ("MK.pl -outfile MK.tmp/$transcript.MKout.tmp1.txt -pol pol -outgroup $outgroup1 -outgroup $outgroup2 -ingroup $ingroup -dir ./MK.tmp/");
					system ("MK_trim.pl MK.tmp/$transcript.MKout.tmp1.txt > MK.tmp/$transcript.MKout.tmp2.txt");	
				
				# run MK.pl with "unpol" option	
				} elsif ($pol eq "unpol") {
					system ("MK.pl -outfile MK.tmp/$transcript.MKout.tmp1.txt -pol unpol -outgroup $outgroup -ingroup $ingroup -dir ./MK.tmp/");
					system ("MK_trim.pl MK.tmp/$transcript.MKout.tmp1.txt > MK.tmp/$transcript.MKout.tmp2.txt");
				}
				# Clean-up seqs from MK tmp directory for the next sequence
				system ("rm MK.tmp/$transcript.fa ");
			}	
			
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
			system ("mv $ctl_dir/* .");
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
	
	
			### 2. BRANCH MODELS 
			system ("sed -i.bak 's/model = 0/model = 2/g' codeml.ctl");	
				# 2a. AMR
				# edit tree
				system ("sed -i.bak 's/amr /amr \#1/g' ANLV.tree");
				# edit codeml.ctl 
				system ("sed -i.bak 's/paml.out/$transcript.br_amr.out/g' codeml.ctl");
				# run & fetch
				system ("codeml");
				my $A_br_lnL = `grep -w -m 1 'lnL' $transcript.br_amr.out | awk '{print \$5}'`;
				# back to normal
				system ("sed -i.bak 's/amr \#1/amr /g' ANLV.tree");
				system ("sed -i.bak 's/$transcript.br_amr.out/paml.out/g' codeml.ctl");
		
				# 2b. NOV
				# edit tree
				system ("sed -i.bak 's/nov /nov \#1/g' ANLV.tree");
				# edit codeml.ctl 
				system ("sed -i.bak 's/paml.out/$transcript.br_nov.out/g' codeml.ctl");
				# run & fetch
				system ("codeml");
				my $N_br_lnL = `grep -w -m 1 'lnL' $transcript.br_nov.out | awk '{print \$5}'`;
				# back to normal
				system ("sed -i.bak 's/nov \#1/nov /g' ANLV.tree");
				system ("sed -i.bak 's/$transcript.br_nov.out/paml.out/g' codeml.ctl");
				
				# 2c. LUM
				# edit tree
				system ("sed -i.bak 's/lum /lum \#1/g' ANLV.tree");
				# edit codeml.ctl 
				system ("sed -i.bak 's/paml.out/$transcript.br_lum.out/g' codeml.ctl");
				# run & fetch
				system ("codeml");
				my $L_br_lnL = `grep -w -m 1 'lnL' $transcript.br_lum.out | awk '{print \$5}'`;
				# back to normal
				system ("sed -i.bak 's/lum \#1/lum /g' ANLV.tree");
				system ("sed -i.bak 's/$transcript.br_lum.out/paml.out/g' codeml.ctl");
		
				# 2d. VIRq
				# edit tree
				system ("sed -i.bak 's/vir /vir \#1/g' ANLV.tree");
				# edit codeml.ctl 
				system ("sed -i.bak 's/paml.out/$transcript.br_vir.out/g' codeml.ctl");
				# run & fetch
				system ("codeml");
				my $V_br_lnL = `grep -w -m 1 'lnL' $transcript.br_vir.out | awk '{print \$5}'`;
				# back to normal
				system ("sed -i.bak 's/vir \#1/vir /g' ANLV.tree");
				system ("sed -i.bak 's/$transcript.br_vir.out/paml.out/g' codeml.ctl");
		
				# 2e. AMR-NOV
				# edit tree
				system ("sed -i.bak 's/nov ),/nov ) \$1,/g' ANLV.tree");
				# edit codeml.ctl 
				system ("sed -i.bak 's/paml.out/$transcript.br_amr-nov.out/g' codeml.ctl");
				# run & fetch
				system ("codeml");
				my $AN_br_lnL = `grep -w -m 1 'lnL' $transcript.br_amr-nov.out | awk '{print \$5}'`;
				# back to normal
				system ("sed -i.bak 's/nov ) \$1,/nov ),/g' ANLV.tree");
				system ("sed -i.bak 's/$transcript.br_amr-nov.out/paml.out/g' codeml.ctl");
				
				# 2f. AMR-NOV-LUM
				# edit tree
				system ("sed -i.bak 's/lum ),/lum ) \$1,/g' ANLV.tree");
				# edit codeml.ctl 
				system ("sed -i.bak 's/paml.out/$transcript.br_amr-nov-lum.out/g' codeml.ctl");
				# run & fetch
				system ("codeml");
				my $ANL_br_lnL = `grep -w -m 1 'lnL' $transcript.br_amr-nov-lum.out | awk '{print \$5}'`;
				# back to normal
				system ("sed -i.bak 's/lum ) \$1,/lum ),/g' ANLV.tree");
				system ("sed -i.bak 's/$transcript.br_amr-nov-lum.out/paml.out/g' codeml.ctl");
			system ("sed -i.bak 's/model = 2/model = 0/g' codeml.ctl");
			system ("mv $transcript.br_*.out PAML.output");
	
	
			### 3. SITE MODELS
			system ("sed -i.bak 's/NSsites = 0/NSsites = 1 2 7 8/g' codeml.ctl");
				# edit codeml.ctl
				system ("sed -i.bak 's/paml.out/$transcript.site_1278.out/g' codeml.ctl");
				# run & fetch
				system ("codeml");
				my $M1_lnL = `grep 'lnL' $transcript.site_1278.out | head -2 | tail -1 | awk '{print \$5}' | sed 's/)://g' `;
				my $M2_lnL = `grep 'lnL' $transcript.site_1278.out | head -3 | tail -1 | awk '{print \$5}' | sed 's/)://g' `;
				my $M7_lnL = `grep 'lnL' $transcript.site_1278.out | head -4 | tail -1 | awk '{print \$5}' | sed 's/)://g' `;
				my $M8_lnL = `grep 'lnL' $transcript.site_1278.out | head -5 | tail -1 | awk '{print \$5}' | sed 's/)://g' `;
				#back to normal
				system ("sed -i.bak 's/$transcript.site_1278.out/paml.out/g' codeml.ctl");
			system ("sed -i.bak 's/NSsites = 1 2 7 8/NSsites = 0/g' codeml.ctl");
			system ("mv $transcript.site_1278.out PAML.output");
		
			### 4. BRANCH-SITE MODELS
			system ("sed -i.bak 's/NSsites = 0/NSsites = 2/g' codeml.ctl");
			system ("sed -i.bak 's/model = 0/model = 2/g' codeml.ctl");
			system ("sed -i.bak 's/fix_omega = 0/fix_omega = 1/g' codeml.ctl");
				# 4a. AMR
				# edit tree
				system ("sed -i.bak 's/amr /amr \#1/g' ANLV.tree");
				# edit codeml.ctl 
				system ("sed -i.bak 's/paml.out/$transcript.brSt_amr.H0.out/g' codeml.ctl");
				### run and get H0 info
				system ("codeml");
				my $A_brSt_H0 = `grep -w -m 1 'lnL' $transcript.brSt_amr.H0.out | awk '{print \$5}'`;
				### swith to H1
				system ("sed -i.bak 's/$transcript.brSt_amr.H0.out/$transcript.brSt_amr.H1.out/g' codeml.ctl");
				system ("sed -i.bak 's/fix_omega = 1/fix_omega = 0/g' codeml.ctl");
				### run and get H1 info
				system ("codeml");
				my $A_brSt_H1 = `grep -w -m 1 'lnL' $transcript.brSt_amr.H1.out | awk '{print \$5}'`;
				### back top normal
				system ("sed -i.bak 's/$transcript.brSt_amr.H1.out/paml.out/g' codeml.ctl");
				system ("sed -i.bak 's/fix_omega = 0/fix_omega = 1/g' codeml.ctl");
				system ("sed -i.bak 's/amr \#1/amr /g' ANLV.tree");

				# 4b. NOV
				# edit tree
				system ("sed -i.bak 's/nov /nov \#1/g' ANLV.tree");
				# edit codeml.ctl 
				system ("sed -i.bak 's/paml.out/$transcript.brSt_nov.H0.out/g' codeml.ctl");
				### run and get H0 info
				system ("codeml");
				my $N_brSt_H0 = `grep -w -m 1 'lnL' $transcript.brSt_nov.H0.out | awk '{print \$5}'`;
				### swith to H1
				system ("sed -i.bak 's/$transcript.brSt_nov.H0.out/$transcript.brSt_nov.H1.out/g' codeml.ctl");
				system ("sed -i.bak 's/fix_omega = 1/fix_omega = 0/g' codeml.ctl");
				### run and get H1 info
				system ("codeml");
				my $N_brSt_H1 = `grep -w -m 1 'lnL' $transcript.brSt_nov.H1.out | awk '{print \$5}'`;
				### back top normal
				system ("sed -i.bak 's/$transcript.brSt_nov.H1.out/paml.out/g' codeml.ctl");
				system ("sed -i.bak 's/fix_omega = 0/fix_omega = 1/g' codeml.ctl");
				system ("sed -i.bak 's/nov \#1/nov /g' ANLV.tree");
				
				# 4c. LUM
				# edit tree
				system ("sed -i.bak 's/lum /lum \#1/g' ANLV.tree");
				# edit codeml.ctl 
				system ("sed -i.bak 's/paml.out/$transcript.brSt_lum.H0.out/g' codeml.ctl");
				### run and get H0 info
				system ("codeml");
				my $L_brSt_H0 = `grep -w -m 1 'lnL' $transcript.brSt_lum.H0.out | awk '{print \$5}'`;
				### swith to H1
				system ("sed -i.bak 's/$transcript.brSt_lum.H0.out/$transcript.brSt_lum.H1.out/g' codeml.ctl");
				system ("sed -i.bak 's/fix_omega = 1/fix_omega = 0/g' codeml.ctl");
				### run and get H1 info
				system ("codeml");
				my $L_brSt_H1 = `grep -w -m 1 'lnL' $transcript.brSt_lum.H1.out | awk '{print \$5}'`;
				### back top normal
				system ("sed -i.bak 's/$transcript.brSt_lum.H1.out/paml.out/g' codeml.ctl");
				system ("sed -i.bak 's/fix_omega = 0/fix_omega = 1/g' codeml.ctl");
				system ("sed -i.bak 's/lum \#1/lum /g' ANLV.tree");
		
				# 4d. VIR
				# edit tree
				system ("sed -i.bak 's/vir /vir \#1/g' ANLV.tree");
				# edit codeml.ctl 
				system ("sed -i.bak 's/paml.out/$transcript.brSt_vir.H0.out/g' codeml.ctl");
				### run and get H0 info
				system ("codeml");
				my $V_brSt_H0 = `grep -w -m 1 'lnL' $transcript.brSt_vir.H0.out | awk '{print \$5}'`;
				### swith to H1
				system ("sed -i.bak 's/$transcript.brSt_vir.H0.out/$transcript.brSt_vir.H1.out/g' codeml.ctl");
				system ("sed -i.bak 's/fix_omega = 1/fix_omega = 0/g' codeml.ctl");
				### run and get H1 info
				system ("codeml");
				my $V_brSt_H1 = `grep -w -m 1 'lnL' $transcript.brSt_vir.H1.out | awk '{print \$5}'`;
				### back top normal
				system ("sed -i.bak 's/$transcript.brSt_vir.H1.out/paml.out/g' codeml.ctl");
				system ("sed -i.bak 's/fix_omega = 0/fix_omega = 1/g' codeml.ctl");
				system ("sed -i.bak 's/vir \#1/vir /g' ANLV.tree");
		
				# 4e. AMR-NOV
				# edit tree
				system ("sed -i.bak 's/nov ),/nov ) \$1,/g' ANLV.tree");
				# edit codeml.ctl 
				system ("sed -i.bak 's/paml.out/$transcript.brSt_amr-nov.H0.out/g' codeml.ctl");
				### run and get H0 info
				system ("codeml");
				my $AN_brSt_H0 = `grep -w -m 1 'lnL' $transcript.brSt_amr-nov.H0.out | awk '{print \$5}'`;
				### swith to H1
				system ("sed -i.bak 's/$transcript.brSt_amr-nov.H0.out/$transcript.brSt_amr-nov.H1.out/g' codeml.ctl");
				system ("sed -i.bak 's/fix_omega = 1/fix_omega = 0/g' codeml.ctl");
				### run and get H1 info
				system ("codeml");
				my $AN_brSt_H1 = `grep -w -m 1 'lnL' $transcript.brSt_amr-nov.H1.out | awk '{print \$5}'`;
				### back top normal
				system ("sed -i.bak 's/$transcript.brSt_amr-nov.H1.out/paml.out/g' codeml.ctl");
				system ("sed -i.bak 's/fix_omega = 0/fix_omega = 1/g' codeml.ctl");
				system ("sed -i.bak 's/nov ) \$1,/nov ),/g' ANLV.tree");
				
				# 4e. AMR-NOV-LUM
				# edit tree
				system ("sed -i.bak 's/lum ),/lum ) \$1,/g' ANLV.tree");
				# edit codeml.ctl 
				system ("sed -i.bak 's/paml.out/$transcript.brSt_amr-nov-lum.H0.out/g' codeml.ctl");
				### run and get H0 info
				system ("codeml");
				my $ANL_brSt_H0 = `grep -w -m 1 'lnL' $transcript.brSt_amr-nov-lum.H0.out | awk '{print \$5}'`;
				### swith to H1
				system ("sed -i.bak 's/$transcript.brSt_amr-nov-lum.H0.out/$transcript.brSt_amr-nov-lum.H1.out/g' codeml.ctl");
				system ("sed -i.bak 's/fix_omega = 1/fix_omega = 0/g' codeml.ctl");
				### run and get H1 info
				system ("codeml");
				my $ANL_brSt_H1 = `grep -w -m 1 'lnL' $transcript.brSt_amr-nov-lum.H1.out | awk '{print \$5}'`;
				### back top normal
				system ("sed -i.bak 's/$transcript.brSt_amr-nov-lum.H1.out/paml.out/g' codeml.ctl");
				system ("sed -i.bak 's/lum ) \$1,/lum ),/g' ANLV.tree");
			system ("sed -i.bak 's/NSsites = 2/NSsites = 0/g' codeml.ctl");
			system ("sed -i.bak 's/model = 2/model = 0/g' codeml.ctl");
			system ("mv $transcript.brSt_*.out PAML.output");
			system ("sed -i.bak 's/$transcript.phy/paml.phy/g' codeml.ctl");	
			
			
			################## COMBINE OUTPUTS ######################
			my $outfile = "$transcript.comp_paml.tmp.txt";
			my $out;
			open($out, ">$outfile") or die "can't open $outfile: $!\n";
			print $out 			"$transcript\t$H0_omega\t$H0_lnL\t$A_br_lnL\t$N_br_lnL\t$V_br_lnL\t$AN_br_lnL\t$M1_lnL\t$M2_lnL\t$M7_lnL\t$M8_lnL\t$A_brSt_H0\t$A_brSt_H1\t$N_brSt_H0\t$N_brSt_H1\t$V_brSt_H0\t$V_brSt_H1\t$AN_brSt_H0\t$AN_brSt_H1\n";
		

 			system ("rm $transcript.*P.afa");
 			system ("rm $transcript.phy");
			system ("cat $transcript.comp_paml.tmp.txt | awk -v OFS='\t' '\$1=\$1' | tr '\n' '\t' | perl -p -e 's/\t\$/\n/g' > PAML.output/$transcript.comp_paml.txt");
			system ("rm $transcript.comp_paml.tmp.txt");
			system ("mv ANLV.tree codeml.ctl $ctl_dir");
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
		if ($mk_test_flag) {
			system ("echo 'TRANSCRIPT\tNS_POLY\tS_POLY\tNS_FIX\tS_FIX\tcodons\tfinal_NI\talpha\tfinal_FET' | cat - MK.tmp/*.MKout.tmp2.txt > $output_dir/MKout.txt");
			system ("rm -r MK.tmp");
		}
		
		if ($calculate_KaKs_flag) {
			system ("cat kaks.tmp/*KaKs.txt | sed '/TRANSCRIP/d' > kaks.tmp/All.tmp.KaKs.txt");
			system ("echo 'TRANSCRIPT\tSEQ1\tSEQ2\tKa\tKs\tKa/Ks\tPROT_PERCENTID\tCDNA_PERCENTID' | cat - kaks.tmp/All.tmp.KaKs.txt > $output_dir/FINAL.KaKs.txt");
			system ("rm -r KaKs.tmp");
		}
		
		if ($run_PAML_flag){
				system ("echo 'TRANSCRIPT\tpaml_dNdS\tH0_lnL\tA_br_lnL\tN_br_lnL\tV_br_lnL\tAN_br_lnL\tM1_lnL\tM2_lnL\tM7_lnL\tM8_lnL\tA_brSt_H0\tA_brSt_H1\tN_brSt_H0\tN_brSt_H1\tV_brSt_H0\tV_brSt_H1\tAN_brSt_H0\tAN_brSt_H1' | cat - PAML.output/*.comp_paml.txt > PAML.output/comp_paml.all.txt");
				system ("sed -i.bak '/^\$/d' PAML.output/comp_paml.all.txt");
				system ("rm *bak");
				system ("mv PAML.output $output_dir/");
				system ("rm r* 2N* 4* lnf");
		}	
}


