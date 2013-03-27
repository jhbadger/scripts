#!/usr/bin/env perl 

use strict;
use DBI;

@ARGV == 2 || die ("Usage: $0 mysql-server database\n");

my $dbh = DBI->connect("dbi:mysql:host=$ARGV[0];database=$ARGV[1]",
    "access","access") or die "Connect Error:",  $DBI::errstr;

my $sql = "SELECT kingdom, phylum, class, ord, family, genus, species FROM classification";

my $request = $dbh->prepare($sql) or die "Prepare Error: ", $dbh->errstr;
$request->execute or die "Execute Error: ", $dbh->errstr;

my %counts = ();

while (my @data = $request->fetchrow_array()) {
  my $taxon = "";
  for(my $i=0; $i< @data; $i++) {
    if ($data[$i] ne "Mixed" && $data[$i] ne "Undefined") {
      $taxon = $data[$i];
    }
  }
  if ($taxon ne "") {
    $counts{$taxon} += 1;
  }
}

foreach my $key (sort {$counts{$b} <=> $counts{$a} } keys %counts) {
    printf("%s\t%d\n", $key, $counts{$key});
}


