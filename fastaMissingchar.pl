#!/usr/bin/perl

my $usage="\nUsage: $0 [-h] [-m char] [fastaFileName1 ...]\n".
    "  -h: help\n".
    "  -m: missing character\n".
    "Print out the name of sequences with characters other than ATGC-.\n".
    "If -m is specified, the ambiguous characters are repleced with the\n".
    "specified character.  e.g. -m '?' will place ? to the ambigous characters.\n" .
    "If multiple files are given, sequences in all files are marged.  If no \n".
    "argument is given, it will take STDIN as the input\n";

our($opt_h, $opt_m);

use Bio::SeqIO;

use Getopt::Std;
getopts('hm:') || die "$usage\n";
die "$usage\n" if (defined($opt_h));

my $format = "fasta";
my @seqArr = ();

@ARGV = ('-') unless @ARGV;
while (my $file = shift) {
    my $seqio_obj = Bio::SeqIO->new(-file => $file, -format => $format);
    while (my $seq = $seqio_obj->next_seq()) {
	push(@seqArr, $seq);
    }
}

#@seqArr = sort { $a->id() cmp $b->id() } @seqArr;

foreach my $s (@seqArr) {
    my $thisSeq = $s->seq();
    my $ambig = AmbiguousChar($thisSeq);
    if ($ambig ne "") {
	print STDERR $s->id(), "\t$ambig\n";
	if (defined($opt_m)) {
	    $thisSeq = ReplaceAmbiguousChar($thisSeq, $opt_m);
	    $s->seq($thisSeq);
	}
    }
}

if (defined($opt_m)) {
    my $seqOut = Bio::SeqIO->new(-fs => \*STDOUT, -format => $format);
    foreach my $s (@seqArr) {
	$seqOut->write_seq($s);
    }
}
exit;


sub AmbiguousChar {
    my $string = shift;
    $string =~ s/[ATGCatgc]//g;

    $string =~ s/\s+//g;
    return $string;
}

sub ReplaceAmbiguousChar {
    my ($string, $char) = @_;
    $string =~ s/[^ATGCatgc]/$char/g;
    return $string;
}

