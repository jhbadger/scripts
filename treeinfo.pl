#!/usr/bin/env perl
#
# newick2treeinfo.pl 
# written by John P. McCrow
# 
# Revision history: 
#

use strict;

my $intree = shift;
my $ininfo = shift;

if($intree) {
} else {
  die "Usage: $0 [Newick Tree file] ([Taxa Info file])\n";
}

my $treestr;
my $maxnodeid=0;
my $maxtaxid=0;
my %nodechild;
my %nodelen;
my %nodeboot;
my %nodesubtaxa;
my %nctype;
my %maxchild;
my %taxname;
my %taxlen;
my $node;
my $tax;
my $boot;
my ($id, $name, $assign);
my %taxtext;
my %taxassign;
my $i;
my $haslengths=0; #Are lengths given

open(IN, $intree) or die "Unable to open file $intree\n";
while(<IN>) {
  chomp;
  $treestr .= $_;
}
close(IN);

if($ininfo) {
  open(IN, $ininfo) or die "Unable to open file $ininfo\n";
  while(<IN>) {
    chomp;
    ($id, $name, $assign) = split(/\t/);
    $taxtext{$id} = $name;
    $taxassign{$id} = $assign;
  }
  close(IN);
}

if($treestr =~ /^\(([^;]+)\);$/) {
  $nodelen{0} = 0;
  parsetreestr(1, length($1), 0, 1);
} else {
  die "Unexpected tree format, should be (...);\n";
}

#If no lengths are given, output unit lengths
print ">nodes\n";
foreach $node (sort {$a<=>$b} keys %maxchild) {
  for($i=1;$i<=$maxchild{$node};$i++) {
    print join("\t", ($node, $haslengths?$nodelen{$node}:1, 0+$nodeboot{$node}, $nodesubtaxa{$node}, $nctype{$node.".".$i}, $nodechild{$node.".".$i}))."\n";
  }
}
print ">taxa\n";
foreach $tax (sort {$a<=>$b} keys %taxlen) {
  if(!defined($taxtext{$taxname{$tax}})) {
    $taxtext{$taxname{$tax}} = $taxname{$tax};
  }
  print join("\t", ($tax, $haslengths?$taxlen{$tax}:1, "r", $taxtext{$taxname{$tax}}, $taxassign{$taxname{$tax}}))."\n";
}

while(<DATA>) {
  print $_;
}

###

sub parsetreestr($$$$) {
  my $strbegin = shift;
  my $strend = shift;
  my $parentnode = shift;
  my $child = shift;
  my $i;
  my $j;
  my $ci;
  my $cj;
  my $nodeid;
  my $len;
  my $name;
  my $parencount;
  my $readtoken = 1;
  my $readlength = 0;
  my $subtaxa = 0;

  for($i=$strbegin; $i<=$strend; $i++) {
    $ci = substr($treestr, $i, 1);

    if($ci eq "(") {

      #find closing paren
      $parencount=1;
      for($j=$i+1; $parencount>0 && $j<=$strend; $j++) {
	$cj = substr($treestr, $j, 1);
	if($cj eq "(") {
	  $parencount++;
	} elsif($cj eq ")") {
	  $parencount--;
	}
      }
      if($j-2 > $i+1 && $j-2 < $strend && $parencount == 0) { #recursively evaluate subtree
	$maxnodeid++;

	$nodeid = $maxnodeid;
	$subtaxa += parsetreestr($i+1, $j-2, $maxnodeid, 1);
	
	$i = $j-1;

	$nctype{$parentnode.".".$child} = "n";

	$readlength = 0;
	$readtoken = 0;
      
      } else {
	die "No matching closing parenth for open parenth at position $i\n";
      }

    } elsif($ci eq ",") {

      if($nctype{$parentnode.".".$child} eq "t") {
	if(length($len)>0) {
	  $haslengths=1;
	}
	if(length($name)>0) {
	  $maxtaxid++;
	  $taxname{$maxtaxid} = $name;
	  $nodechild{$parentnode.".".$child} = $maxtaxid;
	  $taxlen{$maxtaxid} = $len;
	  $subtaxa++;
	} else {
	  die "Expected token and length of child $child under node $parentnode\n";
	}
      } elsif($nctype{$parentnode.".".$child} eq "n") {
	if(length($nodeid)>0) {
	  $nodechild{$parentnode.".".$child} = $nodeid;
	  $nodelen{$nodeid} = $len;
	  $nodeboot{$nodeid} = $boot;
	} else {
	  die "Expected subtree and length of child $child under node $parentnode\n";
	}
      } else {
	die "Unexpected ',' in node $parentnode\n";
      }

      $readtoken = 1;
      $readlength = 0;
      $name = "";
      $len = "";
      $boot = "";
      $child++;

    } elsif($ci eq ":") {

      $readlength = 1;
      $readtoken = 0;

    } elsif($readlength && $ci =~ /[\.\d\-e]/) {

      $len .= $ci;

    } elsif($readtoken) {

      $name .= $ci;

      $nctype{$parentnode.".".$child} = "t";

    } elsif(!$readtoken && !$readlength && $ci =~ /[\.\d\-]/) {
      
      $boot .= $ci;

    } else {
      die "Unexpected character \"$ci\" at position $i\n";
    }

  }

  if(length($len)>0) {
    $haslengths=1;
  }
  if($nctype{$parentnode.".".$child} eq "t") {
    if(length($name)>0) {
      $maxtaxid++;
      $taxname{$maxtaxid} = $name;
      $nodechild{$parentnode.".".$child} = $maxtaxid;
      $taxlen{$maxtaxid} = $len;
      $subtaxa++;
    } else {
      die "Expected token and length of child $child under node $parentnode\n";
    }
  } elsif($nctype{$parentnode.".".$child} eq "n") {
    if(length($nodeid)>0) {
      $nodechild{$parentnode.".".$child} = $nodeid;
      $nodelen{$nodeid} = $len;
      $nodeboot{$nodeid} = $boot;
    } else {
      die "Expected subtree and length of child $child under node $parentnode\n";
    }
  } else {
    die "Unexpected end of node $parentnode\n";
  }

  $maxchild{$parentnode} = $child;
  $nodesubtaxa{$parentnode} = $subtaxa;

  return $subtaxa;
}

__DATA__
>annotations
plot	radial	0.1	0	1	#000000
rad	5	0
font	1.5	Verdana	#000000
labs	1	2	#000000
leg	0.1	0.9	0.9

# All annotations are tab-delimited with type first, followed by parameters
#	offsets are measured in fraction of the default or the total viewable area
#	taxonomies must be semi-colon delimited
# plot (tree type, x-offset, y-offset, tree zoom, tree color)
# rad  (open circle deg, radial tree rotation deg)
# font (font scale, font family, font color)
# labs (label offset, label scale, label color)
# leg  (legend size/text, legend x-offset, legend y-offset)

# rootn (node#) plots the subtree below node#
#ex. rootn	10

# roott (taxonomy) plots the subtree that includes all of the taxonomy given
#ex. roott	Bacteria; Proteobacteria; Gammaproteobacteria

# ttclr (tax type, font color, highlight color) color of specified taxon type
#ex. ttclr	q	#FF0000	#FFFF00

#htax (taxonomy, color) highlights taxa matching taxonomy with color
#dx. htax	bacteria	#DDFFDD

#hst (node#, color) highlights subtree below node# with color
#ex. hst	64	#AAFFAA

#label (tax#, text, color) places label of text near tax#, color is optional
#ex. label	793	Archaea
#ex. label	770	Eukaryota	#FF0000

#tvchart (scale, grid baselines (0,1), color1, color2, second chart position(top,side) ) plots taxval values in 1 or 2 histograms
#ex. tvchart	1	1	#0000FF	#00FF00	top

#tvgrid (value, chart(1,2) ) plots a gray line at value according to chart 1 or 2, default is 1
#ex. tvgrid	10	2

#taxval (taxid, value1, value2) values to plot for each taxon, second value is optional
#ex. taxval	123.45
#ex. taxval	123.45	67.89

#absize (abundance radius scale) factor to scale abundance radius values
#ex. absize	10

#abund (node#, value, color) draws a filled in circle at node# with radius relative to value, scaled by largest abund value
