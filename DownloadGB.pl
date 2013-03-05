#!/usr/bin/env perl

use LWP::Simple;
$url="http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&rettype=fasta&retmode=text&id=";
$list = "";
$count = 0;
while(<>) {
	$_.chomp;
	$list .= $_.",";
	$count +=1;
	if ($count > 100) {
		$count = 0;
		getstore($url.$list, "tmp.txt");
		open(T, "tmp.txt");
		while (<T>) {
			print $_;
		}
		$list = "";
		sleep 0.1
	}
}
chop($list);
if ($count > 0) {
		$count = 0;
		getstore($url.$list, "tmp.txt");
		open(T, "tmp.txt");
		while (<T>) {
			print $_;
		}
	}
