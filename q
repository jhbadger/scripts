#!/usr/bin/env perl
use strict;

my $Iam = $ENV{"USER"};
my $option = shift;

my $qstatcmd;
my $qstr;
my ($pid, $pr, $jname);
my ($usr, $stat, $date, $time, $queue);
my $qtype;
my %jnumberedjobs;
my %collapsejobs;
my $jbase;

my %qtypecount;
my %usercount;
my %statcount;
my %queuestat;
my %usrstat;
my $maxusrlen;
my @allstatsordered;

my %myjobs;
my %mystats;
my %myjobstats;
my $mymaxjnamelen;
my @mystatsordered;

if($option) {
  $qstatcmd = "qstat -u \"*\"";
} else {
  $qstatcmd = "qstat -u ".$Iam;
}

print localtime().": ".$qstatcmd."\n";

$qstr = `$qstatcmd`;
  
foreach (split(/\n/,$qstr)) {
  if (/^\s*(\d+)\s([\d\.]+)\s(..........)\s+/) {
    ($pid, $pr, $jname) = ($1, $2, $3);
    ($usr, $stat, $date, $time, $queue) = split(/\s+/,$');    #');
    if($queue =~ /^[a-z]/i) {
    } else {
      $queue = "queued";
    }
    ($qtype) = split(/\./, $queue);
    $qtypecount{$qtype}++;
    $usercount{$usr}++;
    $statcount{$stat}++;
    $queuestat{$qtype.".".$stat}++;
    $usrstat{$usr.".".$stat}++;
    if(length($usr) > $maxusrlen) {
      $maxusrlen = length($usr);
    }

    if($usr eq $Iam) {
	if($jname =~ /\d+\s*$/) {
	    $jbase = "".$`."*";    #`;
	    $jnumberedjobs{$jbase}++;
	    $myjobstats{$jbase.".".$stat}++;
	}
	$myjobs{$jname}++;
	$mystats{$stat}++;
	$myjobstats{$jname.".".$stat}++;
	if (length($jname) > $mymaxjnamelen) {
	    $mymaxjnamelen = length($jname);
	}
    }
    
  }
}

if(scalar(keys %myjobs)>0) {
    foreach $jname (sort keys %myjobs) {
	if($jname =~ /\d+\s*$/) {
	    $jbase = "".$`."*";    #`;
	    if($jnumberedjobs{$jbase}>1) {
		$collapsejobs{$jbase}=$jnumberedjobs{$jbase};
		$collapsejobs{$jname}=0;
	    }
	}
    }

    foreach $jname (keys %collapsejobs) {
	if($collapsejobs{$jname}>0) {
	    $myjobs{$jname}=$collapsejobs{$jname};
	} else {
	    delete $myjobs{$jname};
	}
    }
}

push(@allstatsordered, "r");
push(@allstatsordered, "qw");
foreach $stat (sort keys %statcount) {
  if ($stat eq "r" || $stat eq "qw") {
  } else {
    push(@allstatsordered, $stat);
  }
}

push(@mystatsordered, "r");
push(@mystatsordered, "qw");
foreach $stat (sort keys %mystats) {
  if ($stat eq "r" || $stat eq "qw") {
  } else {
    push(@mystatsordered, $stat);
  }
}


if($option && scalar(keys %qtypecount)>0) {
  print "### All queues ~ users\n";
  printf "%".$maxusrlen."s\tTotal", "  ";
  foreach $stat (@allstatsordered) {
    print "\t".$stat;
  }
  print "\n";

  foreach $qtype (sort keys %qtypecount) {
    printf "%".$maxusrlen."s\t%d", $qtype, $qtypecount{$qtype};
    foreach $stat (@allstatsordered) {
      print "\t".(0+$queuestat{$qtype.".".$stat});
    }
    print "\n";
  }

  print "~\n";

  foreach $usr (sort keys %usercount) {
    printf "%".$maxusrlen."s\t%d", $usr, $usercount{$usr};
    foreach $stat (@allstatsordered) {
      print "\t".(0+$usrstat{$usr.".".$stat});
    }
    print "\n";
  }
  print "\n";
}

if(scalar(keys %myjobs)>0) {
  if($option) {
    print "### User ".$Iam."\n";
  }
  printf "%".$mymaxjnamelen."s\tTotal", "  ";
  foreach $stat (@mystatsordered) {
    print "\t$stat";
  }
  print "\n";
  
  foreach $jname (sort keys %myjobs) {
    printf "%".$mymaxjnamelen."s\t%d", $jname, $myjobs{$jname};
    foreach $stat (@mystatsordered) {
      print "\t".(0+$myjobstats{$jname.".".$stat});
    }
    print "\n";
  }
}
