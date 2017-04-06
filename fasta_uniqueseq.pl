#!/usr/bin/perl

use strict;
use warnings;

BEGIN 
{
    unshift @INC, '/Users/yazahmed/Programs/FASTA_Unique_Sequences_1.0/FASTA_Unique_Sequences_1.0/modules';
}

use Jls::Software;
use Jls::MultiFasta;
use Jls::Fasta;
use Jls::File;
use Jls::Search;

my $input;
my $output;
my $unique = 1;

my $str;

if ($#ARGV == -1) { Jls::Software::error ('', -1); }

for (my $arg = 0; $arg <= $#ARGV; $arg++)
{
    my $s = substr ($ARGV [$arg], 0, 2);

    if ($s eq '-r') { $unique = 0; }
    elsif ($arg == $#ARGV) { Jls::Software::error ($ARGV [$arg], 1); }
    elsif ($s eq '-i')
    {
        if ($input) { Jls::Software::error ($ARGV [$arg], 2); }
        
        $input = $ARGV [++$arg];
        $str = Jls::File::file_to_string ($input);
    }
    elsif ($s eq '-o')
    {
        if ($output) { Jls::Software::error ($ARGV [$arg], 2); }
        $output = $ARGV [++$arg];
        open (STDOUT, ">$output") or 
            die "Redirection of STDOUT to $output failed.";
    }
    else { Jls::Software::error ($ARGV [$arg], 0); }
}

if (! $input) { die 'No input file was specified with the -i argument.'; }

my $in_fastas;
my $out_fastas;

$in_fastas = Jls::Fasta::string2fastas ($str); 

if ($unique) # unique
{
    $out_fastas = Jls::MultiFasta::fastas2multi_fastas ($in_fastas);
}
else # de-unique
{
    my @fas = ();
    
    foreach my $fasta (@$in_fastas)
    {
        my @deflines = split ('>', substr ($fasta->{'def'}, 1));
        
        foreach my $defline (@deflines)
        {
            my $record = '>' . $defline . "\n" . $fasta->{'seq'};
            push (@fas, new Jls::Fasta ($record))
        }
    }
    
    my @a = sort { $a->{'def'} cmp $b->{'def'} } @fas;

    $out_fastas = \@a;
}

$str = Jls::Fasta::fastas2string ($out_fastas);
print STDOUT $str;
close STDOUT;

sub name
{
    return 'fasta_uniqueseq.pl 1.0';
}

sub args
{
my $msg = '
The program uniques the sequences in a FASTA file.
It concatenates deflines that share sequences. 
It can also reverse the concatenation.
See the URL below about criteria pertaining to shared sequences.

Program options are arguments with \'-\' followed by a letter.
An option requiring further input(s) appears with square brackets.
The default for an option, if any, is indicated in parentheses.
An option with no input always has a default; it lacks square brackets.

-i [input FASTA file] 
-o [output alphabetized FASTA file, for overwriting] (default : STDOUT)
-r  reverse the concatenation

Please see the URL
http://www.ncbi.nlm.nih.gov/CBBresearch/Spouge/html.ncbi/fasta/list.html
for further information.';

    return $msg;
}
