#!/usr/bin/perl

my $usage="Usage: $0 [-h] [inputFile]\n".
    " Read in fasta file, and convert lowercase base to 'n'.\n" .
    " This may be useful to remove unreliable sites of contigs from consed.\n".
    " STDIN is used as the input if no fastaFile is given\n";


my $sep = "\t";  # if you use tab in the sequence name, change this to
                 # other characters such as ","

use Getopt::Std;
getopts('hf:') || die "$usage\n";
die "$usage\n" if (defined($opt_h));

@ARGV = ('-') unless @ARGV; # take STDIN when no arg.

## read in seq data
my $seqFile = shift @ARGV;

my @dat = ReadInFASTAnoCaseChange($seqFile);

foreach my $line (@dat) {
    my ($name, $seq) = split /$sep/, $line;
    print ">$name\n";
    $seq =~ s/[atgc]/N/g;
    print "$seq\n";
}

exit (0);

sub ReadInFASTAnoCaseChange {
    my $infile = shift;
    my @line;
    my $i = -1;
    my @result = ();
    my @seqName = ();
    my @seqDat = ();

    open (INFILE, "<$infile") || die "Can't open $infile\n";

    while (<INFILE>) {
        chomp;
        if (/^>/) {  # name line in fasta format
            $i++;
            s/^>\s*//; s/^\s+//; s/\s+$//;
            $seqName[$i] = $_;
            $seqDat[$i] = "";
        } else {
            s/^\s+//; s/\s+$//;
            s/\s+//g;                  # get rid of any spaces
            next if (/^$/);            # skip empty line
            s/u/t/g;                  # change U to T
            s/U/T/g;                  # change U to T
            $seqDat[$i] = $seqDat[$i] . $_;
        }

        # checking no occurence of internal separator $sep.
        die ("ERROR: \"$sep\" is an internal separator.  Line $. of " .
             "the input FASTA file contains this charcter. Make sure this " .
             "separator character is not used in your data file or modify " .
             "variable \$sep in this script to some other character.\n")
            if (/$sep/);

    }
    close(INFILE);

    foreach my $i (0..$#seqName) {
        $result[$i] = $seqName[$i] . $sep . $seqDat[$i];
    }
    return (@result);
}
