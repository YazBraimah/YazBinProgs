#!/usr/bin/perl

=head1 SYNOPSIS

use to get log likelihoods:

paml_TMP.pl <gene_list> <blast_list_file> <cds_file> <codeml_ctrl>

  
=cut

use strict;
use warnings;
use Pod::Usage;
use Bio::Tools::Run::Phylo::PAML::Codeml;
use Bio::Tools::Run::Alignment::Clustalw;
 
# for projecting alignments from protein to R/DNA space
use Bio::Align::Utilities qw(aa_to_dna_aln);
# for input of the sequence data
use Bio::SeqIO;
use Bio::AlignIO;
 
my $input_file = $ARGV[0];
my $blast_file = $ARGV[1];
my $CDS_file = $ARGV[2];
my $codeml_ctrl_file = $ARGV[3];
my ($help);

my $aln_factory = Bio::Tools::Run::Alignment::Clustalw->new;

pod2usage(1) if($help);
pod2usage("No files given!\n")  if ((!$input_file || !$blast_file));

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
system ("grep $int $blast_file | tr '\t' '\n' | sed '/$int/d' > $int.list");
system ("faSomeRecords $CDS_file $int.list $int.fa");

system ("sed -i.bak '/amr/ s/>/>amr /' $int.fa");
system ("sed -i.bak '/lum/ s/>/>lum /' $int.fa");
system ("sed -i.bak '/nov/ s/>/>nov /' $int.fa");
system ("sed -i.bak '/vir/ s/>/>vir /' $int.fa");
system ("sed -i.bak '/moj/ s/>/>moj /' $int.fa");
system ("sed -i.bak '/gri/ s/>/>gri /' $int.fa");
system ("rm *bak");

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
          warn("provided a CDS sequence with a stop codon, PAML will choke!");
          exit(0);
    }
    # Tcoffee can't handle '*' even if it is trailing
    $pseq =~ s/\*//g;
    $protein->seq($pseq);
    push @prots, $protein;
}
 
if( @prots < 2 ) {
    warn("Need at least 2 CDS sequences to proceed");
    exit(0);
}

open(OUT, ">$int.KaKs.txt") ||  die("cannot open output align_output for writing");

# Align the sequences with clustalw
my $aa_aln = $aln_factory->align(\@prots);
# project the protein alignment back to CDS coordinates
my $dna_aln = aa_to_dna_aln($aa_aln, \%seqs); 
my @each = $dna_aln->each_seq();

########################################################################################################
###################### THIS SEGMENT TO WRITE OUT SEQS AND RUN PAML #####################################

my $out_dna = Bio::AlignIO->new(-file => ">$int.aln_dna.tmp.afa" ,
							-format => 'fasta');
$out_dna -> write_aln($dna_aln);

my $out_aa = Bio::AlignIO->new(-file => ">$int.aln_aa.tmp.afa" ,
							-format => 'fasta');
$out_aa -> write_aln($aa_aln);

system ("fastaSortByName.pl $int.aln_dna.tmp.afa > $int.aln_dna.afa");
system ("fastaSortByName.pl $int.aln_aa.tmp.afa > $int.aln_aa.afa");
system ("sed -i.bak '/>/ s/\\/.*//g' $int.*.afa");
system ("paml_prep.pl $int.aln_aa.afa $int.aln_dna.afa -nogap -output paml > $int.phy");
system ("perl -pi -e 's/paml.phy/$int.phy/g' $codeml_ctrl_file");
system ("perl -pi -e 's/paml.out/$int.out/g' $codeml_ctrl_file");
system ("codeml $codeml_ctrl_file");
system ("perl -pi -e 's/$int.phy/paml.phy/g' $codeml_ctrl_file");
system ("perl -pi -e 's/$int.out/paml.out/g' $codeml_ctrl_file");

########################################################################################################
########################################################################################################
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
system ("cat *KaKs.txt | sed '/TRANSCRIP/d' > All.tmp.KaKs.txt");
system ("echo 'TRANSCRIPT\tSEQ1\tSEQ2\tKa\tKs\tKa/Ks\tPROT_PERCENTID\tCDNA_PERCENTID' | cat - All.tmp.KaKs.txt > FINAL.KaKs.txt");
system ("mkdir final_results");
system ("mv FINAL.KaKs.txt final_results");
system ("rm *list *.fa *txt *bak");

close FILE1;