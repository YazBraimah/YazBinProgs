#!/usr/bin/perl

=head1 SYNOPSIS

 EXON-ext.pl <bam_file> <gff_file> <genome> <orientation_file> <species_name> <isoform_name> <gene_name> <XLOC-ID>

	bam_file: obvious!
	
	gff_file: annotation file with following tab separated values:gene_name/coordinates/orientation/"exon"/gene-name_exon-number.

	genome: obvious!

	orientation file: two-column file with gene name and orientation (from gff_file: awk '{print $1, "\t", $3}' gff_file)

	species_name[STR]: three-letter species extenstion

	isoform_name: name of isoform as it appears on first columns of GFF file.
	
	gene_name: GJ ID or Cufflinks gene name (CUFF###)
  
=cut


use strict;
use warnings;
use Bio::SeqIO;
use Pod::Usage;

#Define command line argument variables:
my $bam_in = $ARGV[0];
my $gff_file = $ARGV[1];
my $genome = $ARGV[2];
my $tcons = $ARGV[5];
my $species = $ARGV[4];
my $gene_id = $ARGV[6];
my $xloc_id = $ARGV[7];
my ($help);


pod2usage(1) if($help);
pod2usage("No files given!\n")  if ((!$bam_in || !$gff_file));


#open GFF file and process:
open (FILE2, $gff_file) || die  ("Cannot open file \n");
while ($gff_file = <FILE2>){
chomp $gff_file;

my @gene_bits = split(/\s+/, $gff_file);
my $new_gene = $gene_bits[0];
my $coordinates = $gene_bits[1];
my $type = $gene_bits[3];
my $exon_num = $gene_bits[4];

my @chase = ();
push(@chase, $new_gene);

foreach my $arm (@chase){

	if ($arm eq $tcons ) {


		system ("mkdir $gene_id\_$xloc_id\_$tcons/");

        system ("samtools mpileup -Iq20 -r $coordinates -gf $genome $bam_in > ./$gene_id\_$xloc_id\_$tcons/$tcons.bcf");

		if($arm eq $tcons && $type eq 'exon') {

             system ("bcftools view -g ./$gene_id\_$xloc_id\_$tcons/$tcons.bcf  - | vcfutils.pl vcf2fq -Q20 | seqtk seq -A - | seqtk cutN - | sed 's/>/>$coordinates/g' > $gene_id\_$xloc_id\_$tcons/$coordinates\_$exon_num.fasta");

            }

		}
	}
}
system ("cat ./$gene_id\_$xloc_id\_$tcons/*.fasta > ./$gene_id\_$xloc_id\_$tcons/$tcons.fa");
system ("convert.pl ./$gene_id\_$xloc_id\_$tcons/$tcons.fa > ./$gene_id\_$xloc_id\_$tcons/$tcons\_masked.fa");

system ("rm ./$gene_id\_$xloc_id\_$tcons/*.fasta");

system ("perl -pi -e 's/.txt//g' ./$gene_id\_$xloc_id\_$tcons/$tcons\_masked.fa");





	

close(FILE2);


#########  PROCESS THE FASTA FILE  ###############

my $filename = $ARGV[1];


### gff file as the same file in the first part of the script ###

open (FILE4, $filename) || die ("Could not get file \n");
while ($filename = <FILE4>){
chomp $filename;

my @name_bits = split(/\s+/, $filename);
my $final_gene = $name_bits[0];
my $new_pos = $name_bits[4];
my $mol_type = $name_bits[3];


my @names = ();
push (@names, $final_gene);

my @coord = ();
push (@coord, $new_pos);

		foreach my $new_gene(@names){

	if ($new_gene eq $tcons){						



foreach my $add (@coord){



					
					

my $fasta_file = "./$gene_id\_$xloc_id\_$tcons/$tcons\_masked.fa";  # ID to extract

my $sample_final = "./$gene_id\_$xloc_id\_$tcons/$tcons\_masked.CD";



local $/ = "\n>";  # read by FASTA record

open (FILE6, $fasta_file) || die ("Could not open file \n");

open (FINAL, ">$sample_final") or die ("Could not open file \n");


my @fasta = ();

while ($fasta_file = <FILE6>) {
    chomp $fasta_file;

        $fasta_file =~ /^>*(\S+)/;  # parse ID as first word in FASTA header
        $fasta_file =~ s/^>*.+\n//;  # remove FASTA header
        $fasta_file =~ s/\n//g;  # remove endlines

push (@fasta, $fasta_file);


}

print FINAL ">$gene_id\_$xloc_id\_$tcons\_$species\n@fasta\n";




		}


	}
	
	}

}	

close (FILE4);
close (FILE6);
close(FILE8);
close(FINAL);


### process the files for final alignment ###

my $genes = $ARGV[3];


open (FILE9, $genes) || die ("Could not get file \n");
while ($genes = <FILE9>){
chomp $genes;

my @prev_bits = split(/\s+/, $genes);
my $prev_gene = $prev_bits[0];
my $strand = $prev_bits[1];


my @prev = ();
push (@prev, $prev_gene);

        foreach my $new_gene(@prev){

	if ($new_gene eq $tcons && $strand eq '+'){


system ("cat ./$gene_id\_$xloc_id\_$tcons/*.CD > ./$gene_id\_$xloc_id\_$tcons/$tcons\_align_exons.fasta");

system ("fastaMissingchar.pl -m N ./$gene_id\_$xloc_id\_$tcons/$tcons\_align_exons.fasta > ./$gene_id\_$xloc_id\_$tcons/$tcons\_exons.FASTA");



					}


    elsif ($new_gene eq $tcons && $strand eq '-'){

system ("cat ./$gene_id\_$xloc_id\_$tcons/*.CD > ./$gene_id\_$xloc_id\_$tcons/$tcons\_align_exons.fasta");

system ("fastaMissingchar.pl -m N ./$gene_id\_$xloc_id\_$tcons/$tcons\_align_exons.fasta > ./$gene_id\_$xloc_id\_$tcons/$tcons\_exons.FASTA");




my $seqin = Bio::SeqIO->new(-file => "./$gene_id\_$xloc_id\_$tcons/$tcons\_exons.FASTA",
                            -format => "fasta");
my $seqout = Bio::SeqIO->new(-file => ">./$gene_id\_$xloc_id\_$tcons/$tcons\_exons_rev.FASTA",
                            -format => "fasta");



while(my $seq = $seqin->next_seq){

my $new_seq = $seq->revcom;

$seqout->write_seq($new_seq);

		}



		}

	}

}

system ("rm ./$gene_id\_$xloc_id\_$tcons/*.CD");
system ("rm ./$gene_id\_$xloc_id\_$tcons/*.fa");
system ("rm ./$gene_id\_$xloc_id\_$tcons/*.fasta*");
system ("rm ./$gene_id\_$xloc_id\_$tcons/*.bcf");
###############################
system ("mkdir Gene_exons");
system ("mv $gene_id\_$xloc_id\_$tcons ./Gene_exons");





close (FILE8);
