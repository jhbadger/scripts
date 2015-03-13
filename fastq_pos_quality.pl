#!/usr/bin/env perl
use strict;

my @fqfiles = @ARGV;

unless(scalar(@fqfiles) > 0) {
    die "Usage: $0 [FASTQ file(s)]\n";
}

my $maxp = 0;
my %posave;

foreach my $infq (@fqfiles) {
    
    my %poscount = ();
    my %possum = ();

    if($infq =~ /\.gz$/) {
	open(IN, "gunzip -c $infq 2>/dev/null |");
    } else {
	open(IN, $infq) or die "Unable to open file $infq\n";
    }

    print STDERR "Reading file: $infq ...\n";

    my $rec = 0;
    my $sn = 0;
    do {
	$sn++;
	if(defined($_ = <IN>)) {
	    chomp;
	    if($sn == 4) {
		my @qlist = split(//, $_);
		for(my $i=0; $i<scalar(@qlist); $i++) {
		    my $v = ord($qlist[$i])-33;
		    $possum{$i} += ($v * 0.01);
		    $poscount{$i}++;
		}
	    }
	}
	if($sn >= 4) {
	    $sn = 0;
	    $rec++;
	}
    } until(eof IN || $rec >= 1000000);

    foreach my $p (keys %poscount) {
	$posave{$infq}{$p} = 100.0 * ($possum{$p} / $poscount{$p});
	if($p > $maxp) {
	    $maxp = $p;
	}
    }
}

print "position";
foreach my $f (@fqfiles) {
    print "\t".$f;
}
print "\n";

for(my $i=0; $i<=$maxp; $i++) {
    printf "%d", $i+1;
    foreach my $f (@fqfiles) {    
	printf "\t%.2f", $posave{$f}{$i};
    }
    print "\n";
}
