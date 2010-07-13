#!/usr/bin/env perl -w

use strict;
use DBI;

my $start = time;

my $infile = $ARGV[0];
my $dbh = DBI->connect('dbi:mysql:host=mysql-lan-pro;database=phylodb','amoustafa','changeme') or die "Connect Error:",  $DBI::errstr;

my $i = 0;
while (my $line = <>) {
    $i++;
    chomp($line);

    my @array = split(/\t/, $line);

    my $query = $array[0];
    my $phylodb_seq_id = $array[1];

    print $query, "\t", $phylodb_seq_id;

    my $sqly = qq { select ifnull(t2.taxon_id, ''), ifnull(t2.species,'') from proteins t1 inner join contigs t2 on t1.contig_name = t2.name where t1.name = '$phylodb_seq_id' } ;
    my $sthy = $dbh->prepare($sqly) or die "Prepare Error: ", $dbh->errstr;
    $sthy->execute or die "Execute Error: ", $dbh->errstr;
    my $taxon_id = '';
    my $species = '';
    $sthy->bind_columns(\$taxon_id, \$species);
    $sthy->fetch();
    $sthy->finish();

    $species =~ s/\s/_/g;
    $species =~ s/\.//g;

    print "\t", $taxon_id, "\t", $species;


    my $sqlx = qq { select t1.hit, t2.definition from blast t1 inner join swiss t2 on t1.hit = t2.swiss where t1.source = 'phylodb' and t1.target = 'swiss' and t1.blast = 'blastp' and t1.query = '$phylodb_seq_id' };
    my $sthx = $dbh->prepare($sqlx) or die "Prepare Error: ", $dbh->errstr;
    $sthx->execute or die "Execute Error: ", $dbh->errstr;
    my $swiss_hit = '';
    my $swiss_definition = '';
    $sthx->bind_columns(\$swiss_hit, \$swiss_definition);
    $sthx->fetch();
    $sthx->finish();
    print "\t", $swiss_hit, "\t", $swiss_definition;

    
    my $sql0 = qq { select t1.hit, t2.definition from blast t1 inner join kegg t2 on t1.hit = t2.kegg_seq_id where t1.source = 'phylodb' and t1.target = 'kegg' and t1.blast = 'blastp' and t1.query = '$phylodb_seq_id' };
    my $sth0 = $dbh->prepare($sql0) or die "Prepare Error: ", $dbh->errstr;
    $sth0->execute or die "Execute Error: ", $dbh->errstr;
    my $_hit = '';
    my $_definition = '';
    $sth0->bind_columns(\$_hit, \$_definition); 
    $sth0->fetch();
    $sth0->finish();
    print "\t", $_hit, "\t", $_definition;

    
    
    my $sql1 = qq { select ifnull(gos_cluster_id, ''), ifnull(phytax_cluster_id, ''), ifnull(plant,''), ifnull(animal,''), ifnull(ko,''), ifnull(pfam, ''), ifnull(tigrfams, '') from map where phylodb_seq_id = '$phylodb_seq_id' } ;
    my $sth1 = $dbh->prepare($sql1) or die "Prepare Error: ", $dbh->errstr;
    $sth1->execute or die "Execute Error: ", $dbh->errstr;
    my $gos = '';
    my $phytax = '';
    my $plant = '';
    my $animal = '';
    my $kegg = '';
    my $pfam = '';
    my $tigrfams = '';
    $sth1->bind_columns(\$gos, \$phytax, \$plant, \$animal, \$kegg, \$pfam, \$tigrfams);
    $sth1->fetch();
    $sth1->finish();

    print "\t", $gos, "\t", $phytax;

    if ($kegg && $kegg ne '') {

	my $sql2 = qq { select ifnull(t2.ko,''), ifnull(t2.pathway, '') from kegg_entries t1 inner join pathways t2 on t1.ko = t2.ko where t1.entry = '$kegg' } ;
	my $sth2 = $dbh->prepare($sql2) or die "Prepare Error: ", $dbh->errstr;
	$sth2->execute or die "Execute Error: ", $dbh->errstr;
	my $ko = '';
	my $pathway = '';
	$sth2->bind_columns(\$ko, \$pathway);
	$sth2->fetch();
	$sth2->finish();
	
	print "\t", $ko, "\t", $pathway;
    } else {
	print "\t", '' , "\t", '';
    }

    print "\t", $pfam;

    if ($pfam && $pfam ne '') {

        my $sql3 = qq { select ifnull(definition,'') from pfam where pfam = '$pfam' } ;
        my $sth3 = $dbh->prepare($sql3) or die "Prepare Error: ", $dbh->errstr;
        $sth3->execute or die "Execute Error: ", $dbh->errstr;
        my $pfam_definition = '';
        $sth3->bind_columns(\$pfam_definition);
        $sth3->fetch();
        $sth3->finish();

        print "\t", $pfam_definition;

	my $_pfam = $pfam;
	if ($pfam =~ m/(.*)\.\d+/) {
	    $_pfam = $1;
	}

        my $sql3x = qq { select go from pfam2go where pfam = '$_pfam' } ;
        my $sth3x = $dbh->prepare($sql3x) or die "Prepare Error: ", $dbh->errstr;
        $sth3x->execute or die "Execute Error: ", $dbh->errstr;
        my $go = '';
	my $pfam2go = '';
	my $flag = 1;
        $sth3x->bind_columns(\$go);
        while ($sth3x->fetch()) {
	    if ($flag) {
		$pfam2go = $go;
		$flag = 0;
	    } else {
		$pfam2go = $pfam2go . "," . $go;
	    }
	}
        $sth3x->finish();


	print "\t", trim($pfam2go);

    } else {
	print "\t", "", "\t", "";
    }
    
    print "\t", $tigrfams;

    if ($tigrfams && $tigrfams ne '') {

        my $sql4 = qq { select ifnull(definition,'') from tigrfams where tigrfams = '$tigrfams' } ;
        my $sth4 = $dbh->prepare($sql4) or die "Prepare Error: ", $dbh->errstr;
        $sth4->execute or die "Execute Error: ", $dbh->errstr;
        my $tigrfams_definition = '';
        $sth4->bind_columns(\$tigrfams_definition);
        $sth4->fetch();
        $sth4->finish();

        print "\t", $tigrfams_definition;

        my $sql4x = qq { select go from tigrfams2go where tigrfams = '$tigrfams' } ;
        my $sth4x = $dbh->prepare($sql4x) or die "Prepare Error: ", $dbh->errstr;
        $sth4x->execute or die "Execute Error: ", $dbh->errstr;
        my $go = '';
        my $tigrfams2go = '';
        my $flag = 1;
        $sth4x->bind_columns(\$go);
        while ($sth4x->fetch()) {
            if ($flag) {
                $tigrfams2go = $go;
                $flag = 0;
            } else {
                $tigrfams2go = $tigrfams2go . "," . $go;
            }
        }
        $sth4x->finish();

        print "\t", trim($tigrfams2go);

    } else {
        print "\t", "", "\t", "";
    }

    print "\t" , $plant, "\t", $animal;
    print "\n";
}

$dbh->disconnect;

my $end = time;

print STDERR "Elapsed time: " ;
print STDERR $end - $start;
print STDERR " seconds\n";


sub trim {
    my $string = shift;
    $string =~ s/^\s+//;
    $string =~ s/\s+$//;
    return $string;
}
