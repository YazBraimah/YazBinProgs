#!/usr/bin/perl
# get CONSENSUS seq from a BAM file for all transcripts in a transcript.gtf file from cufflinks
# ARGS:
# 1. BAM
# 2. Transcripts.gtf file from cufflinks


=head1 NAME

 get_consensus_bam_batch

=head1 SYNOPSIS

 get_consensus_bam_batch.pl -b file.BAM -t transcripts.GTF

 Options:
 -b - bam file from Tophat
 -t - transcript.gtf file produced by cufflinks
 -h - prints this msg and exits
 
 Examples:
 
 - commons usage:
  get_consensus_bam_batch.pl -b file.BAM -t transcripts.GTF

		
=head1 DESCRIPTION

 This program extracts consensus sequence from a bam file. It requires 
 a modified version of samtools-0.1.18 or higher to generate a pileup file 
 which is parsed for every transcript in the  transcript.gtf produced 
 by cufflinks. That modified version is not skipping alignments which are not
 primary or duplicates. Creates a full pile up for any given region and extracts
 the consensus for that region.

=head1 AUTHOR

 Dimitar Kenanov, <dimitark@aitbiotech.com>

=cut

use strict;
use warnings;
use Getopt::Long;
use List::Util qw(max);
use Pod::Usage;

my ($opts,$help,$man,$bam_in,$trans_in,$trans_out,$pile_out,$err_file,@chroms);

$opts=GetOptions("b=s"=>\$bam_in,
				 "t=s"=>\$trans_in,
				 "h"=>\$help) or pod2usage(2);

pod2usage(1) if($help);
pod2usage("$0: No files given!\n")  if ((!$bam_in || !$trans_in));

#if(!$bam_in || !$trans_in){
#	print "\n Usage:\tget_consensus_bam_batch -b bam -t transcripts\n\n";
#	exit;
#}

$trans_out="$trans_in.v3.fa";
$pile_out="$bam_in.v3.pileup";

# file to show transcript with no start or end
$err_file="$bam_in.log";

my ($fhi,$fho,$fhe,$line,$region,$gid,$rid,$fpkm,$nreads,$found_seqs,$ntrans,@tmp);


open($fho,'>',$trans_out) or die "problems $trans_out: $!\n";
open($fhe,'>',$err_file) or die "problems $err_file: $!\n";

$|=1;

$found_seqs=0;
$ntrans=0;

	open($fhi,'<',$trans_in) or die "problems $trans_in: $!\n";
	while($line=<$fhi>){
		chomp($line);
		if($line=~/\ttranscript\t/){
			$line=~/FPKM\s\"(\S+)\"/;
			$fpkm=$1;
#		print "FPKM:$fpkm:\n";
			if($fpkm > 0){
#			print "FPKM:$fpkm:\n";
				$ntrans++;
				@tmp=split/\t/,$line;
				if($tmp[8]=~/^gene_id\s"(\S+)";/){
					$gid=$1;
#					print "GID1:$gid:$1\n";
				}else{
					$gid='';
				}
			
				if($tmp[8]=~/transcript_id\s\"(\S+)\"/){
					$rid=$1;
				}else{
					$rid='';
				}
			
				if($gid eq ''){				
					$gid=$rid;
#				print "GID2:$gid:\n";
				}elsif($rid eq ''){
					$gid=$gid;
				}else{
					$gid="$gid:$rid";
				}

				my $seq;
				$seq=extract_seq($tmp[0],$tmp[3],$tmp[4]);

			
				if($seq){
					print "\n\tPRINTED FASTA SEQ FOR REGION $tmp[0]:$tmp[3]-$tmp[4]: $gid:";
					print $fho ">$gid; $tmp[0]:$tmp[3]-$tmp[4]\n$seq\n";
					print "THE SEQ:\n$seq\n";
					$found_seqs++;
				}else{
					print $fhe "NO SEQ FOUND FOR REGION: $tmp[0]:$tmp[3]-$tmp[4]\n";
				}
				print "\n";
#			<STDIN>;
			}	
		}
	}
	close $fhi;


print "VALID TRANSCRIPTS: $ntrans :\n";
print "FOUND SEQS: $found_seqs :\n";


close $fho;
close $fhe;

###################
###### SUBS #######
###################

sub extract_seq{
	my ($chr,$start,$end)=@_;
	my ($region,$res,$ret_seq,$n,$i);
	my (@seq,@pileup);
#orig	$region="$chr:$start..$end";
	$region="$chr:$start-$end";
	print "$bam_in;$region;$pile_out\n";
#orig	$res=`bamtools convert -format pileup -in $bam_in -region $region > $pile_out`;
	$res=`samtools_ndsm mpileup -ABQ0 -d 1000000 -r $region $bam_in> $pile_out`;	
	put_in_memory($pile_out,\@pileup,$chr);
#	print "ARR:@{$pileup[0]} \n";
	foreach $i (@pileup){
		process_line(\@$i,\@seq,\$n);
	}
	$ret_seq=join('',@seq);
	return $ret_seq;
}

sub process_line{
	my ($inarr,$arr,$n)=@_;
	my @tmp;
	my ($sa,$sc,$sg,$st);
	my (@a,@c,@g,@t);
	my %svals;
	@a=$$inarr[4] =~/(a)/ig;
	@c=$$inarr[4] =~/(c)/ig;
	@g=$$inarr[4] =~/(g)/ig;
	@t=$$inarr[4] =~/(t)/ig;
	$sa=scalar(@a);
	$sc=scalar(@c);
	$sg=scalar(@g);
	$st=scalar(@t);
	$svals{'a'}=$sa;
	$svals{'c'}=$sc;
	$svals{'g'}=$sg;
	$svals{'t'}=$st;
	
	my $max=0;
	my $letter='';
	my $k;
	foreach $k (keys %svals){
		if($svals{$k}>$max){
			$max=$svals{$k};
			$letter=$k;
		}
	}		
	
#print "MAX: is $letter: $max\n";
	if($max && $letter){
#		$nr=scalar(@tmp);
#		print "TMP:@tmp :$nr\n";
		push @$arr,uc($letter);
	}else{
		$$n++;
	}
}

sub put_in_memory{
	my ($file,$arr,$chr)=@_;
	my ($fh,$line,@tmp);
	open($fh,'<',$file) or die "problems $file:$!\n";
	while($line=<$fh>){
		chomp($line);
		@tmp=split/\t/,$line;
		if($tmp[0] eq $chr){
			push @$arr,[@tmp[0..4]];
		}
	}
	close $fh;
}


