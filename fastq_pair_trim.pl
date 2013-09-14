#!/usr/bin/env perl
use strict;
use Getopt::Long;

my $paired = 0;
my $minq = 0;
my $minlen = 0;
my $outfile;

GetOptions ("q=i" => \$minq,
            "l=i" => \$minlen,
	    "o=s" => \$outfile);

my $help = <<HELP;
FASTQ Quality Trimmer - jmccrow 05/29/2013
Assumes format Phred+33.  Give 2 files to interlace paired data.

Usage: $0 (options) [FASTQ file(s) ...]
    -q int  : Minimum quality score (default: 0)

    -l file : Minimum sequence length (default: 0)

    -o file : Output fasta file (default: STDOUT)

HELP

my $infile1 = shift;
my $infile2 = shift;

unless(length($infile1) > 0) {
    die $help;
}

if(length($infile2) > 0) {
    $paired = 1;
}

if(length($outfile) > 0) {
    open(OUT, ">".$outfile) or die "Unable to write to file $outfile\n";
} else {
    open(OUT, ">&=STDOUT") or die "Unable to write to STDOUT\n";
}

if (index($infile1, ".gz")) {
    open(IN1, "zcat $infile1 |") or die "Unable to open file $infile1\n";
}
else {
  open(IN1, $infile1) or die "Unable to open file $infile1\n";
}


if($paired) {
    if (index($infile2, ".gz")) {
	open(IN2, "zcat $infile2 |") or die "Unable to open file $infile2\n";
    }
    else {
	open(IN2, $infile2) or die "Unable to open file $infile2\n";
    }
}

my $sn = 0;
my %sfline;
do {
    $sn++;
    if(defined($_ = <IN1>)) {
	chomp;
	$sfline{1}{$sn} = $_;
    }
    if($paired && defined($_ = <IN2>)) {
	chomp;
	$sfline{2}{$sn} = $_;
    }
    if($sn >= 4) {
	my $p1 = conv($sfline{1}{1}, $sfline{1}{2}, $sfline{1}{3}, $sfline{1}{4});
	if($paired) { 
	    my $p2 = conv($sfline{2}{1}, $sfline{2}{2}, $sfline{2}{3}, $sfline{2}{4});
	    if(length($p1) > 0 && length($p2) > 0) {
		print OUT $p1.$p2;
	    }
	} else {
	    if(length($p1) > 0) {
		print OUT $p1;
	    }
	}

	%sfline = ();
	$sn = 0;
    }
} until(eof IN1 && (!$paired || eof IN2));


###

sub conv($$$$) {
    my ($l1, $l2, $l3, $l4) = @_;

    my $head = "";
    if($l1 =~ /^\@(.+)$/) {
	$head = $1;
    } else {
	die "Not a valid FASTQ format 1\n";
    }

    my $seq = $l2;

    unless($l3 =~ /^\+/) {
	die "Not a valid FASTQ format 2\n";
    }

    my $qual = $l4;
    unless(length($seq) == length($qual)) {
	die "Not a valid FASTQ format 3\n";
    }

    return qtrim($head, $seq, $qual);
}

sub qtrim($$$) {
    my ($head, $seq, $qual) = @_;
    my $minpos = -1;
    my $maxpos = -1;
    my $trimseq = "";

    my @qlist = split(//, $qual);
    for(my $i=0; $i<scalar(@qlist); $i++) {
	my $v = ord($qlist[$i])-33;
	if($v >= $minq) {
	    if($minpos < 0) {
		$minpos = $i;
	    }
	    $maxpos = $i;
	}
    }
    if($maxpos >= $minpos && $minpos >= 0) {
	$trimseq = substr($seq, $minpos, ($maxpos-$minpos+1));
    }
    my $retseq = "";
    if(length($trimseq) > $minlen) {
	$retseq = ">".$head."\n".$trimseq."\n";
    }
    return $retseq;
}
