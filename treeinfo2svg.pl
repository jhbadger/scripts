#!/usr/bin/env perl
#
# treeinfo2svg.pl 
# written by John P. McCrow
# 
# Revision history: 
#

use strict;
use Math::Trig;

my $infile = shift;
my $layer;
my %svglayerstr;
my $xmult;
my $ymult;
my $tmult;
my $toffx;
my $firstnode;
my $rootnode;
my $rootdist;
my $maxtextlen;
my $maxtaxval1;
my $maxtaxval2;
my $maxnodeabund;
my $numtax;
my $rcx;
my $rcy;
my $viewwidth;
my $viewheight;
my $maxfont;
my $minfont;
my $maxx;
my $linewidth;
my $fgcolor;
my $fontsize;
my $taxfontsize;
my $labfontsize;
my $titlefontsize;
my $t;
my $v;
my $hmcolstitl;
my $hastaxval1;
my $hastaxval2;
my @hmcolnames;
my %taxvalgridlines1;
my %taxvalgridlines2;
my %labeltaxclr;
my %labeltaxstr;
my %taxslot;
my %slottax;
my %taxcolor;
my %hmtaxvals;
my %nodelen;
my %nodeboot;
my %nodetaxa;
my %nodect;
my %cid;
my %taxlen;
my %taxtype;
my %taxname;
my @leghstrs;
my %leghfg;
my %leghbg;
my %nodeabund;
my %nodeabundclr;
my %taxassign;
my %taxval1;
my %taxval2;
my %ttfgclr;
my %ttbgclr;
my %hmvclr;
my $firstdataline = 1;

#Default Parameters:
my $plotstyle = 2; #1=normal, 2=radial
my $radialspace = 10;
my $treerotation = 0;
my $aligntext = 1;
my $treecolor="#000000";
my $fontcolor="#000000";
my $labelcolor=$fontcolor;
my $fontfam = "Verdana";
my $legsize = "0.1"; #Reference size legend: used as both the text, and line length
my $treezoom = 1;
my $fontzoom = 2;
my $labelfontmult = 1;
my $xoffset = 0;
my $yoffset = 0;
my $labeloffsetmult = 1;
my $legxoffset = 0.1;
my $legyoffset = 0.1;
my $abundmaxradius = 10;
my $titletext = "";
my $titlefontmult = 4;
my $tvchartscale = 1;
my $tvdefaultgrid = 0;
my $tvchartcolor1 = "#0000FF";
my $tvchartcolor2 = "#00FF00";
my $tvchart2pos = 1;
my $mintvlen = 0.005; # Prevents ultra-small tvchart lines, which are often misinterpreted when converting to pdf

my $hstlayer = 1;
my $treelayer = 2;
my $abundlayer = 3;

if($infile) {
} else {
  die "Usage: $0 [Tree Info File]\n";
}

readtreeinfo();
calcscaleparams();

#Add tree and abundance values to layers
#drawsubtree must be run before labels and annotations, because they rely on the placement of taxa to be known
drawsubtree($rootnode, $rootdist, 0);

dlegend();
dtitle();

foreach $t (keys %labeltaxstr) {
  if(length($labeltaxclr{$t})>0) {
    dlabel($taxslot{$t}, $labeltaxstr{$t}, $labeltaxclr{$t});
  } else {
    dlabel($taxslot{$t}, $labeltaxstr{$t}, $labelcolor);
  }
}

if($plotstyle == 1 && scalar(%hmtaxvals)>0) {
  drawheatmap();
}

#Add hightlighted subtrees to layer
$layer = $hstlayer;
drawhighlightbands();

if($maxtaxval1>0) {
  if($tvdefaultgrid) {
    #Draw default grid lines at the base of each taxval ring that has data
    dgridline(0,1);
    if($maxtaxval2>0) {
      dgridline(0,2);
    }
  }
  foreach $v (keys %taxvalgridlines1) {
    dgridline($v,1);
  }
  foreach $v (keys %taxvalgridlines2) {
    dgridline($v,2);
  }
}

#Print SVG in ordered layers
printheader();
foreach $layer (sort {$a<=>$b} keys %svglayerstr) {
  print $svglayerstr{$layer};
}
printfooter();

###
# Subroutines:
# Naming convention for drawing subs is the following: 
#    draw- high level tree walking, or other calculations, before calling d- subs
#    d- lower level drawing of individual objects in tree space with conversion to absolute printing space
#    print- direct printing of absolute space objects in SVG

sub calcscaleparams() {
  #Calculate scale parameters
  $numtax = getnumtaxa($rootnode);
  $viewwidth = 300000;
  $viewheight = 300000;
  $maxfont = $viewheight/100;
  $minfont = $viewheight/1000;
  $maxx = getmaxx($rootnode, $rootdist);
  $linewidth=int($viewwidth/3000);
  if ($plotstyle == 1) {
    $xmult = ($treezoom*$viewwidth*0.95)/($maxx*1.5);
    $ymult = ($viewheight*0.95)/($numtax+2);
    $rcx = ($viewwidth*$xoffset) + ($viewwidth*0.01);
    $rcy = ($viewheight*$yoffset) + ($viewwidth*0.01);
    $taxfontsize = $fontzoom*($ymult/2);
  } else {
    $rootdist = 0;
    $taxfontsize = $fontzoom*(($viewheight/$numtax)<$minfont?$minfont:(($viewheight/$numtax)>$maxfont?$maxfont:$viewheight/$numtax));
    $xmult = ($treezoom*$viewwidth*0.95)/($maxx*3);
    $ymult = (360-$radialspace) / $numtax;
    $rcx = ($viewwidth*$xoffset) + ($viewwidth/2);
    $rcy = ($viewwidth*$yoffset) + ($viewwidth/2);
  }
  $labfontsize = $labelfontmult*$taxfontsize;
  $titlefontsize = $titlefontmult*$taxfontsize;
  $fontsize = $taxfontsize;
  $toffx = $viewwidth/1500;
  $tmult = $fontsize*0.6;

  $layer = $treelayer;
  $fgcolor = $treecolor;
}

sub readtreeinfo() {
  my ($n, $l, $b, $t, $nt, $cn, $m);
  my ($type, $x, $y, $x1, $x2, $zt, $c, $r);
  my ($str, $zf, $f, $s, $c1, $c2, $tt);
  my $namestr;
  my $taxstr;
  my $readtype;
  my @params;
  my @values;

  open(IN, $infile) or die "Unable to open file $infile\n";
  while (<IN>) {
    chomp;
    if (/^>/) {
      ($readtype) = split(/\s/, $');
      $readtype =~ tr/A-Z/a-z/;
    } else {
      if ($readtype eq "nodes") {

	($n, $l, $b, $t, $nt, $cn) = split(/\t/);
	if($t =~ /[tn]/ && length($cn) < 1) { #No bootstrap values given, shift values to compensate
	  $cn = $nt;
	  $nt = $t;
	  $t = $b;
	}

	if (length($firstnode)==0) {
	  $firstnode = $n;
	  $rootnode = $n;
	}
	$nodelen{$n}=$l;
	$nodeboot{$n}=$b;
	$nodetaxa{$n}=$t;
	push(@{$nodect{$n}}, $nt);
	push(@{$cid{$n}}, $cn);

      } elsif ($readtype eq "taxa") {

	($n, $l, $t, $namestr, $taxstr) = split(/\t/);
	$taxlen{$n}=$l;
	$taxtype{$n}=$t;
	$taxname{$n}=$namestr;
	$taxassign{$n}=$taxstr;
	if (length($namestr) > $maxtextlen) {
	  $maxtextlen = length($namestr);
	}

      } elsif ($readtype eq "annotations") {

	($t, @params) = split(/\t/);

	if ($t =~ /^\#/) {
	  # Skip comment lines
	} elsif ($t eq "plot") {
	  ($type, $x, $y, $zt, $c) = @params;
	  if ($type =~ /^r/i) {
	    $plotstyle = 2;
	  } else {
	    $plotstyle = 1;
	  }
	  $xoffset = 0+$x;
	  $yoffset = 0+$y;
	  if (length($zt)>0) {
	    $treezoom = 0+$zt;
	  }
	  if (length($c)>0) {
	    $treecolor = $c;
	  }

	} elsif ($t eq "rad") {
	  ($a, $r) = @params;
	  if (length($a)>0 && $a>=0 && $a<360) {
	    $radialspace = $a;
	  }
	  if (length($r)>0 && $r>=0 && $r<360) {
	    $treerotation = $r;
	  }
	} elsif ($t eq "font") {
	  ($zf, $f, $s, $c) = @params;
	  if (length($zf)>0) {
	    $fontzoom = 0+$zf;
	  }
	  if (length($f)>0) {
	    $fontfam = $f;
	  }
	  if (length($c)>0) {
	    $fontcolor = $c;
	  }

	} elsif ($t eq "leg") {
	  ($str, $x, $y) = @params;
	  if (length($str)>0) {
	    $legsize = $str;
	    if (length($x)>0) {
	      $legxoffset = $x;
	    }
	    if (length($y)>0) {
	      $legyoffset = $y;
	    }	  
	  }

	} elsif ($t eq "legh") {
	  ($str, $c1, $c2) = @params;
	  push(@leghstrs, $str);
	  $leghfg{$str} = $c1;
	  $leghbg{$str} = $c2;
      
	} elsif ($t eq "labs") {
	  ($x1, $x2, $c) = @params;
	  if (length($x1)>0) {
	    $labeloffsetmult = $x1;
	  }
	  if (length($x2)>0) {
	    $labelfontmult = $x2;
	  }
	  if (length($c)>0) {
	    $labelcolor = $c;
	  }

	} elsif ($t eq "absize") { # Max Abundance Size (absize)
	  ($x) = @params;
	  if (length($x)>0) {
	    $abundmaxradius = $x;
	  }

	} elsif ($t eq "htax") { #Highlight Taxonomy (htax)
	  ($str, $c) = @params;
	  $str =~ s/\s*;\s*/;/g;
	  drawtaxhighlight($rootnode, $str, $c);

	} elsif ($t eq "hst") { #Highlight SubTree (hst)
	  ($n, $c, $m) = @params;
	  if ($m eq "-") {
	    colorsubtreeminus($rootnode, $n, $c);
	  } else {
	    colorsubtree($n, $c);
	  }

	} elsif ($t eq "label") { #Label centered on taxon
	  ($n, $str, $c) = @params;
	  $labeltaxstr{$n} = $str;
	  if (length($c)>0) {
	    $labeltaxclr{$n} = $c;
	  }

	} elsif ($t eq "abund") { #Abundance values on nodes
	  ($n, $a, $c) = @params;
	  $nodeabund{$n} = $a;
	  $nodeabundclr{$n} = $c;
	  if ($a > $maxnodeabund) {
	    $maxnodeabund = $a;
	  }
	} elsif ($t eq "taxval") { #Taxon values to plot as bar
	  my ($n, $v1, $v2) = @params;
	  if(length($v1)>0) {
	    $hastaxval1 = 1;
	    $taxval1{$n}=$v1;
	    if ($v1 > $maxtaxval1) {
	      $maxtaxval1 = $v1;
	    }
	  }
	  if(length($v2)>0) {
	    $hastaxval2 = 1;
	    $taxval2{$n}=$v2;
	    if ($v2 > $maxtaxval2) {
	      $maxtaxval2 = $v2;
	    }
	  }

	} elsif ($t eq "tvchart") {
	  my ($s, $g, $c1, $c2, $p) = @params;
	  if($s > 0) {
	    $tvchartscale = $s;
	  }
	  if($g) {
	    $tvdefaultgrid = 1;
	  }
	  if(length($c1) > 0) {
	    $tvchartcolor1 = $c1;
	  }
	  if(length($c2) > 0) {
	    $tvchartcolor2 = $c2;
	  }
	  if($p =~ /side/ || $p == 2) {
	    $tvchart2pos = 2;
	  }

	} elsif ($t eq "tvgrid") {
	  my ($v, $c) = @params;
	  if($c == 2) {
	    $taxvalgridlines2{$v}=1;
	  } else {
	    $taxvalgridlines1{$v}=1;
	  }

	} elsif ($t eq "ttclr") { #Taxon Type Color, fg and bg colors for taxa of given type
	  ($tt, $c1, $c2) = @params;
	  $ttfgclr{$tt} = $c1;
	  $ttbgclr{$tt} = $c2;

	} elsif ($t eq "rootn") { #Display tree only below given node
	  ($n) = @params;
	  $rootnode = $n;

	} elsif ($t eq "roott") { #Display tree only below given node
	  ($a) = @params;
	  $n = findtaxnode(0, $a);
	  if ($n >= 0) {
	    $rootnode = $n;
	  }

	} elsif ($t eq "titl") { #Title of plot
	  ($str, $n) = @params;
	  $titletext = $str;
	  if ($n > 0) {
	    $titlefontmult = $n;
	  }
	} elsif ($t eq "hmdata") { #heatmap data
	  foreach $str (@params) {
	    ($v, $c) = split(/:/, $str);
	    if (length($v)>0 && length($c)>0) {
	      $hmvclr{$v} = $c;
	    }
	  }
	}

      } elsif ($readtype eq "data") {

	($tt, @values) = split(/\t/);
	if ($firstdataline) {
	  $hmcolstitl = $tt;
	  @hmcolnames = @values;
	} else {
	  if (length($tt)>0 && $tt =~ /^\d+$/) {
	    @{$hmtaxvals{$tt}} = @values;
	  }
	}
	$firstdataline = 0;
      }

    }
  }
  close(IN);
}

sub radialx($$) {
  my $x = shift;
  my $y = shift;

  return ($rcx + sin(deg2rad($y*$ymult+$treerotation))*$x);
}

sub radialy($$) {
  my $x = shift;
  my $y = shift;

  return ($rcy + cos(deg2rad($y*$ymult+$treerotation))*$x);
}

sub radialangle($) {
  my $y = shift;

  return (90-($y*$ymult+$treerotation));
}

sub drawhighlightbands() {
  my $i;
  my $lastclr;
  my $min;
  my $max;

  #Draw color bands for successive same colored taxa
  for($i=0; $i<$numtax; $i++) {
    if(length($lastclr) > 0) {
      if($lastclr eq $taxcolor{$slottax{$i}}) {
	$max = $i;
      } else {
	if($max >= $min && length($lastclr) > 0) {
	  dband($min, $max, $lastclr);
	}
	$min = $i;
	$max = $i;
      }
    } else {
      $min = $i;
      $max = $i;
    }
    $lastclr = "".$taxcolor{$slottax{$i}};
  }
  #Draw last band if necessary
  if($max >= $min && length($lastclr) > 0) {
    dband($min, $max, $lastclr);
  }
  
}

sub drawtaxhighlight($$$) {
  my $n = shift;
  my $taxstr = shift;
  my $color = shift;
  my @ccomp;
  my $ncomp;
  my $i;
  my $c;

  $taxstr =~ tr/A-Z/a-z/;

  for($i=0; $i<scalar(@{$nodect{$n}}); $i++) {
    if(@{$nodect{$n}}[$i] eq "n") {
      push(@ccomp, drawtaxhighlight(@{$cid{$n}}[$i], $taxstr, $color));
    } else {
      $c = $taxassign{@{$cid{$n}}[$i]};
      $c =~ s/\s*;\s*/;/g;
      $c =~ tr/A-Z/a-z/;
      push(@ccomp, $c);
    }
    if($i==0) {
      $ncomp = $ccomp[0];
    } else {
      $ncomp = findcommonassignment($ncomp, $ccomp[$i]);
    }
  }
  
  #Color child node only when current node does not fit taxstr, but child does
  if(length($taxstr) > length($ncomp) || !(substr($ncomp, 0, length($taxstr)) eq $taxstr)) {
    for($i=0; $i<scalar(@{$nodect{$n}}); $i++) {
      if(substr($ccomp[$i], 0, length($taxstr)) eq $taxstr) {
	if(@{$nodect{$n}}[$i] eq "n") {
	  colorsubtree(@{$cid{$n}}[$i], $color);
	} else {
	  $taxcolor{@{$cid{$n}}[$i]} = $color;
	}
      }
    }
  }

  return $ncomp;
}

sub findtaxnode($$) {
  my $n = shift;
  my $t = shift;
  my $nc;
  my $bestnc;
  my $nccount=0;
  my $c;
  my $i;

  $t =~ s/\s*;\s*/;/g;
  $t =~ tr/A-Z/a-z/;

  for($i=0; $i<scalar(@{$nodect{$n}}); $i++) {
    if(@{$nodect{$n}}[$i] eq "n") {
      $nc = findtaxnode(@{$cid{$n}}[$i], $t);
    } elsif(@{$nodect{$n}}[$i] eq "t") {
      $c = $taxassign{@{$cid{$n}}[$i]};
      $c =~ s/\s*;\s*/;/g;
      $c =~ tr/A-Z/a-z/;
      if(substr($c, 0, length($t)) eq $t) {
	$nc = $n;
      } else {
	$nc = -1;
      }
    }
    
    if($nc >= 0) {
      $nccount++;
      $bestnc = $nc;
    }
  }

  if($nccount>1) { # If more than 1 child, then this node or an ancestor is the best inclusive node
    return $n;
  } elsif($nccount==1) { # If only 1 child, then the child is the best inclusive node
    return $bestnc;
  }
  
  return -1;
}

sub getnumtaxa($) {
  my $n = shift;
  my $cn;
  my $count=0;
  my $i;

  for($i=0; $i<scalar(@{$nodect{$n}}); $i++) {
    if(@{$nodect{$n}}[$i] eq "n") {
      $cn = getnumtaxa(@{$cid{$n}}[$i]);
    } elsif(@{$nodect{$n}}[$i] eq "t") {
      $cn=1;
    }
    $count += $cn;
  }

  return $count;
}

sub getmaxx($$) {
  my $n = shift;
  my $px = shift;
  my $x = ($px+$nodelen{$n});
  my $m;
  my $maxx=$x;
  my $i;
  
  for($i=0; $i<scalar(@{$nodect{$n}}); $i++) {
    if(@{$nodect{$n}}[$i] eq "n") {
      if(defined(@{$cid{$n}}[$i])) {
	$m = getmaxx(@{$cid{$n}}[$i], $x);
      } else {
	die "Undefined child node ".($i+1)." under node $n\n";
      }
    } else {
      $m = ($x+$taxlen{@{$cid{$n}}[$i]});
    }
    if($m > $maxx) {
      $maxx = $m;
    }
  }

  return $maxx;
}

sub colorsubtree($$) {
  my $n = shift;
  my $c = shift;
  my $i;
  
  for($i=0; $i<scalar(@{$nodect{$n}}); $i++) {
    if(@{$nodect{$n}}[$i] eq "n") {
      colorsubtree(@{$cid{$n}}[$i], $c);
    } else {
      $taxcolor{@{$cid{$n}}[$i]} = $c;
    }
  }
}

sub colorsubtreeminus($$$) {
  my $n = shift;
  my $notn = shift;
  my $c = shift;
  my $i;
  
  for($i=0; $i<scalar(@{$nodect{$n}}); $i++) {
    if(@{$nodect{$n}}[$i] eq "n") {
      if(@{$cid{$n}}[$i] == $notn) {
      } else {
	colorsubtree(@{$cid{$n}}[$i], $c);
      }
    } else {
      $taxcolor{@{$cid{$n}}[$i]} = $c;
    }
  }
  
}

sub drawsubtree($$$) {
  my $n = shift;
  my $px = shift;
  my $ymin = shift;
  my $y;
  my $yfirst;
  my $ylast;
  my $m;
  my $fc;
  my $ly;
  my $x = ($px+$nodelen{$n});
  my $i;

  for($i=0; $i<scalar(@{$nodect{$n}}); $i++) {
    if(@{$nodect{$n}}[$i] eq "n") {
      $y = drawsubtree(@{$cid{$n}}[$i], $x, $ymin);
      $ymin += $nodetaxa{@{$cid{$n}}[$i]};
    } else {
      $y = drawtax(@{$cid{$n}}[$i], $x, $ymin);
      $ymin++;
    }

    if($i==0) {
      $yfirst = $y;
      $ylast = $y;
    } else {
      $ylast = $y;
    }

  }

  $m = ($yfirst+$ylast)/2;
  
  dvline($x, $yfirst, $ylast);
  dhline($m, $px, $x);
  if($nodeabund{$n} > 0) {
    $ly = $layer;
    $fc = $fgcolor;
    $layer = $abundlayer;
    $fgcolor = $nodeabundclr{$n};
    dcirc($x, $m, sqrt($nodeabund{$n}/$maxnodeabund)*$viewwidth*0.001*$abundmaxradius);
    $fgcolor = $fc;
    $layer = $ly;
  }

  return $m;
}

sub drawtax($$$) {
  my $t = shift;
  my $px = shift;
  my $y = shift;
  my $fc = $fgcolor;

  my $x = ($px+$taxlen{$t});
  $taxslot{$t} = $y;
  $slottax{$y} = $t;

  $fgcolor = $treecolor;
  dhline($y, $px, $x); # Line from parent to terminal
  
  if($aligntext) { # Gray line out to text
    $fgcolor = "#CCCCCC";
    dhline($y, $x, $maxx);
  }

  $fgcolor = $treecolor;
  dcirc($x, $y, $linewidth*2); # Terminal node
  
  if($ttbgclr{$taxtype{$t}}) { # Font and Background by taxon type
    $taxcolor{$t} = $ttbgclr{$taxtype{$t}};
  }
  if($ttfgclr{$taxtype{$t}}) {
    $fgcolor = $ttfgclr{$taxtype{$t}};
  } else {
    $fgcolor = $fontcolor;
  }

  if($aligntext) { # Taxon name
    dhtext($maxx, $y, $taxname{$t});
  } else {
    dhtext($x, $y, $taxname{$t});
  }

  if($hastaxval1) {
    if($hastaxval2) {
      if($taxval1{$t}>0 || $taxval2{$t}>0) {
	if($tvchart2pos == 2) {
	  dbar2_nextto($y, $taxval1{$t}, $taxval2{$t}, $tvchartcolor1, $tvchartcolor2);
	} else {
	  dbar2_ontop($y, $taxval1{$t}, $taxval2{$t}, $tvchartcolor1, $tvchartcolor2);
	}
      }
    } else {
      if($taxval1{$t}>0) {
	dbar($y, $taxval1{$t}, $tvchartcolor1);
      }
    }
  }

  $fgcolor = $fc;

  return $y;
}

sub drawheatmap() {
  my $i;
  my $j;
  my $v;

  dhmplot();
  for($i=0; $i<$numtax; $i++) {
    for($j=0; $j < scalar(@{$hmtaxvals{$slottax{$i}}}); $j++) {
      $v = @{$hmtaxvals{$slottax{$i}}}[$j];
      if(length($hmvclr{$v})>0) {
	dhmbox($j, $i, $hmvclr{$v});
      }
    }
  }

}

sub dhmplot() {
  my $j;
  my $c;
  my $clr;
  my $ax;
  my $ay = $rcy-$ymult;

  for($j=0; $j < scalar(@hmcolnames); $j+=2) {
    if($c==1) {
      $c=0;
      $clr="#DDFFDD";
    } else {
      $c=1;
      $clr="#AAFFAA";
    }
    if($j+1 > scalar(@hmcolnames)-1) {
      dhmrulev($j, scalar(@hmcolnames)-1, $clr);
    } else {
      dhmrulev($j, $j+1, $clr);
    }
  }

  $fontsize = $taxfontsize;
  for($j=0; $j < scalar(@hmcolnames); $j++) {
    $ax = $rcx+($maxx*$xmult)+($maxtextlen*$tmult)+(($j+1.5)*$ymult);
    printtextrot($ax, $ay, -90, $hmcolnames[$j]);
  }

}

sub dhmbox($$$) {
  my $x = shift;
  my $y = shift;
  my $clr = shift;
  my $ax1 = $rcx+($maxx*$xmult)+($maxtextlen*$tmult)+(($x+1)*$ymult);
  my $ay1 = $rcy+($y-0.5)*$ymult;
  my $ax2 = $rcx+($maxx*$xmult)+($maxtextlen*$tmult)+(($x+2)*$ymult);
  my $ay2 = $rcy+($y+0.5)*$ymult;
  my $fc = $fgcolor;

  $fgcolor = $clr;
  printbox($ax1,$ay1,$ax2,$ay2);
  $fgcolor = $fc;
}

sub dhmrulev($$) {
  my $col1 = shift;
  my $col2 = shift;
  my $clr = shift;
  my $ax1 = $rcx+($maxx*$xmult)+($maxtextlen*$tmult)+(($col1+1)*$ymult);
  my $ay1 = $rcy+(-0.5)*$ymult;
  my $ax2 = $rcx+($maxx*$xmult)+($maxtextlen*$tmult)+(($col2+2)*$ymult);
  my $ay2 = $rcy+($numtax+0.5)*$ymult;
  my $fc = $fgcolor;

  $fgcolor = $clr;
  printbox($ax1,$ay1,$ax2,$ay2);
  $fgcolor = $fc;
}

sub dvline($$$) {
  my $x = shift;
  my $y1 = shift;
  my $y2 = shift;
  my $ax1;
  my $ay1;
  my $ax2;
  my $ay2;

  if($plotstyle==2) {
    $ax1 = radialx($x*$xmult, $y1);
    $ay1 = radialy($x*$xmult, $y1);
    $ax2 = radialx($x*$xmult, $y2);
    $ay2 = radialy($x*$xmult, $y2);
    printarc($ax1,$ay1,$ax2,$ay2,$x*$xmult, ((($y2-$y1)*$ymult<180)?0:1) );
  } else {
    $ax1 = $rcx+$x*$xmult;
    $ay1 = $rcy+$y1*$ymult;
    $ax2 = $rcx+$x*$xmult;
    $ay2 = $rcy+$y2*$ymult;
    printline($ax1,$ay1,$ax2,$ay2);
  }
}

sub dhline($$$$) {
  my $y = shift;
  my $x1 = shift;
  my $x2 = shift;
  my $ax1;
  my $ay1;
  my $ax2;
  my $ay2;

  if($plotstyle==2) {
    $ax1 = radialx($x1*$xmult, $y);
    $ay1 = radialy($x1*$xmult, $y);
    $ax2 = radialx($x2*$xmult, $y);
    $ay2 = radialy($x2*$xmult, $y);
  } else {
    $ax1 = $rcx+$x1*$xmult;
    $ay1 = $rcy+$y*$ymult;
    $ax2 = $rcx+$x2*$xmult;
    $ay2 = $rcy+$y*$ymult;
  }
  printline($ax1,$ay1,$ax2,$ay2);

}

sub dhtext($$$) {
  my $x = shift;
  my $y = shift;
  my $str = shift;
  my $ax;
  my $ay;

  $fontsize = $taxfontsize;
  if($plotstyle==2) {
    $ax = radialx($x*$xmult+$toffx+$linewidth*4, $y);
    $ay = radialy($x*$xmult+$toffx+$linewidth*4, $y);
    printtextrot($ax, $ay, radialangle($y), $str);

  } else {
    $ax = $rcx+$toffx+$x*$xmult+$linewidth*4;
    $ay = $rcy+$y*$ymult;
    printtext($ax, $ay, $str);
  }
}

sub dcirc($$$) {
  my $x = shift;
  my $y = shift;
  my $w = shift;
  my $ax;
  my $ay;

  if($plotstyle==2) {
    $ax = radialx($x*$xmult+$toffx, $y);
    $ay = radialy($x*$xmult+$toffx, $y);
    printfillcircle($ax, $ay, $w);

  } else {
    $ax = $rcx+$toffx+$x*$xmult;
    $ay = $rcy+$y*$ymult;
    printfillcircle($ax, $ay, $w);
  }
}

sub dtitle() {
  my $ax;
  my $ay;
  
  $fontsize = $titlefontsize;
  
  $ax = $viewwidth/2;
  $ay = 0;
  
  printtext($ax, $ay, $titletext, "center", "bottom");
}

sub dlegend() {
  my $ax;
  my $ay;
  my $str;
  my $i;
  my $maxlegtext;
  my $fc = $fgcolor;
  my $lw = $linewidth;

  $fontsize = $taxfontsize;

  $ax = $viewwidth*$legxoffset;
  $ay = $viewheight*$legyoffset;

  foreach $str (@leghstrs) {
    if(length($str) > $maxlegtext) {
      $maxlegtext = length($str);
    }
  }

  #Draw any highlighting or text color legend items
  $linewidth = $ymult;
  $i=1;
  foreach $str (@leghstrs) {
    if(length($leghbg{$str}) > 0) {
      $fgcolor = $leghbg{$str};
      printline($ax, $ay+($i*$ymult), $ax+($maxlegtext*$tmult), $ay+($i*$ymult));      
    }
    
    if(length($leghfg{$str}) > 0) {
      $fgcolor = $leghfg{$str};
    } else {
      $fgcolor = $fontcolor;
    } 
    printtext($ax, $ay+($i*$ymult), $str, "left", "center");
    $i++;
  }
  $linewidth = $lw;

  #Draw branch length legend
  if($legsize > 0) {
    $i++;
    $fgcolor = $treecolor;
    printline($ax+($maxlegtext*$tmult)/2, $ay+($i*$ymult), $ax+($maxlegtext*$tmult)/2+($legsize*$xmult), $ay+($i*$ymult));
    $fgcolor = $fontcolor;
    printtext($ax+($maxlegtext*$tmult)/2+($legsize*$xmult)/2, $ay+($i*$ymult)-(0.5*$tmult), $legsize, "center", "bottom");
  }
  
  $fgcolor = $fc;
}

sub dband($$$) {
  my $y1 = shift;
  my $y2 = shift;
  my $color = shift;
  my $fc = $fgcolor;
  my $lw = $linewidth;
  my $ax1;
  my $ay1;
  my $ax2;
  my $ay2;

  #Offset to cover text
  $y1 = $y1-0.5; 
  $y2 = $y2+0.5;

  $fgcolor = $color;
  $linewidth = $maxtextlen*$tmult;

  if($plotstyle==2) {
    $ax1 = radialx($maxx*$xmult+($linewidth/2), $y1);
    $ay1 = radialy($maxx*$xmult+($linewidth/2), $y1);
    $ax2 = radialx($maxx*$xmult+($linewidth/2), $y2);
    $ay2 = radialy($maxx*$xmult+($linewidth/2), $y2);
    printarc($ax1,$ay1,$ax2,$ay2,($maxx*$xmult+($linewidth/2)), ((($y2-$y1)*$ymult<180)?0:1) );
  } else {
    $ax1 = $rcx+($maxx*$xmult+($linewidth/2));
    $ay1 = $rcy+$y1*$ymult;
    $ax2 = $rcx+($maxx*$xmult+($linewidth/2));
    $ay2 = $rcy+$y2*$ymult;
    printline($ax1,$ay1,$ax2,$ay2);
  }

  $fgcolor = $fc;
  $linewidth = $lw;
}

sub dbar($$$) {
  my $y = shift;
  my $v = shift;
  my $color = shift;
  my $fc = $fgcolor;
  my $lw = $linewidth;
  my $ax1;
  my $ay1;
  my $ax2;
  my $ay2;
  my $av = $v/$maxtaxval1;
  $fgcolor = $color;

  if($av > 0 && $av < $mintvlen) {
    $av = $mintvlen;
  }

  if($plotstyle==2) {
    $ax1 = radialx($maxx*$xmult+$maxtextlen*$tmult+$toffx, $y);
    $ay1 = radialy($maxx*$xmult+$maxtextlen*$tmult+$toffx, $y);
    $ax2 = radialx($maxx*$xmult+(($maxtextlen*$tmult)*(1+$av))+$toffx, $y);
    $ay2 = radialy($maxx*$xmult+(($maxtextlen*$tmult)*(1+$av))+$toffx, $y);
    $linewidth = $lw*8;
  } else {
    $ax1 = $rcx+($maxx*$xmult+$maxtextlen*$tmult+$toffx);
    $ay1 = $rcy+$y*$ymult;
    $ax2 = $rcx+($maxx*$xmult+(($maxtextlen*$tmult)*(1+$av))+$toffx);
    $ay2 = $rcy+$y*$ymult;
    $linewidth = $lw*4;
  }
  printline($ax1,$ay1,$ax2,$ay2);

  $fgcolor = $fc;
  $linewidth = $lw;
}

sub dbar2_nextto($$$$$) {
  my $y = shift;
  my $v1 = shift;
  my $v2 = shift;
  my $color1 = shift;
  my $color2 = shift;
  my $fc = $fgcolor;
  my $lw = $linewidth;
  my ($ax11, $ay11, $ax12, $ay12);
  my ($ax21, $ay21, $ax22, $ay22);
  my $inneroffset;
  my $innerringx = ($maxx*$xmult)+($maxtextlen*$tmult)+$toffx;
  my $outerringx = ($maxx*$xmult)+($maxtextlen*$tmult*(1+$tvchartscale))+$toffx;
  my $av1 = $v1/$maxtaxval1;
  my $av2 = $v2/$maxtaxval2;

  if($av1 > 0 && $av1 < $mintvlen) {
    $av1 = $mintvlen;
  }
  if($av2 > 0 && $av2 < $mintvlen) {
    $av2 = $mintvlen;
  }

  if($plotstyle==2) {
    $inneroffset = 0.2;
    $linewidth = $lw * 3;
    $ax11 = radialx($maxx*$xmult+$maxtextlen*$tmult+$toffx, ($y-$inneroffset));
    $ay11 = radialy($maxx*$xmult+$maxtextlen*$tmult+$toffx, ($y-$inneroffset));
    $ax12 = radialx($maxx*$xmult+(($maxtextlen*$tmult)*(1+($tvchartscale*$av1)))+$toffx, ($y-$inneroffset));
    $ay12 = radialy($maxx*$xmult+(($maxtextlen*$tmult)*(1+($tvchartscale*$av1)))+$toffx, ($y-$inneroffset));
    $ax21 = radialx($maxx*$xmult+$maxtextlen*$tmult+$toffx, ($y+$inneroffset));
    $ay21 = radialy($maxx*$xmult+$maxtextlen*$tmult+$toffx, ($y+$inneroffset));
    $ax22 = radialx($maxx*$xmult+(($maxtextlen*$tmult)*(1+($tvchartscale*$av2)))+$toffx, ($y+$inneroffset));
    $ay22 = radialy($maxx*$xmult+(($maxtextlen*$tmult)*(1+($tvchartscale*$av2)))+$toffx, ($y+$inneroffset));
  } else {
    $inneroffset = 0.2;
    $linewidth = $lw * 2;
    $ax11 = $rcx+($maxx*$xmult+$maxtextlen*$tmult+$toffx);
    $ay11 = $rcy+($y-$inneroffset)*$ymult;
    $ax12 = $rcx+($maxx*$xmult+(($maxtextlen*$tmult)*(1+($tvchartscale*$av1)))+$toffx);
    $ay12 = $rcy+($y-$inneroffset)*$ymult;
    $ax21 = $rcx+($maxx*$xmult+$maxtextlen*$tmult+$toffx);
    $ay21 = $rcy+($y+$inneroffset)*$ymult;
    $ax22 = $rcx+($maxx*$xmult+(($maxtextlen*$tmult)*(1+($tvchartscale*$av2)))+$toffx);
    $ay22 = $rcy+($y+$inneroffset)*$ymult;   
  }

  $fgcolor = $color1;
  printline($ax11,$ay11,$ax12,$ay12);
  $fgcolor = $color2;
  printline($ax21,$ay21,$ax22,$ay22);

  $fgcolor = $fc;
  $linewidth = $lw;
}

sub dbar2_ontop($$$$$) {
  my $y = shift;
  my $v1 = shift;
  my $v2 = shift;
  my $color1 = shift;
  my $color2 = shift;
  my $fc = $fgcolor;
  my $lw = $linewidth;
  my ($ax11, $ay11, $ax12, $ay12);
  my ($ax21, $ay21, $ax22, $ay22);
  my $innerringx = ($maxx*$xmult)+($maxtextlen*$tmult)+$toffx;
  my $outerringx = ($maxx*$xmult)+($maxtextlen*$tmult*(1+$tvchartscale))+$toffx;
  my $av1 = $v1/$maxtaxval1;
  my $av2 = $v2/$maxtaxval2;

  if($av1 > 0 && $av1 < $mintvlen) {
    $av1 = $mintvlen;
  }
  if($av2 > 0 && $av2 < $mintvlen) {
    $av2 = $mintvlen;
  }

  if($plotstyle==2) {
    $ax11 = radialx($innerringx, $y);
    $ay11 = radialy($innerringx, $y);
    $ax12 = radialx($innerringx + ($maxtextlen*$tmult*$av1*$tvchartscale), $y);
    $ay12 = radialy($innerringx + ($maxtextlen*$tmult*$av1*$tvchartscale), $y);
    $ax21 = radialx($outerringx, $y);
    $ay21 = radialy($outerringx, $y);
    $ax22 = radialx($outerringx + ($maxtextlen*$tmult*$av2*$tvchartscale), $y);
    $ay22 = radialy($outerringx + ($maxtextlen*$tmult*$av2*$tvchartscale), $y);
    $linewidth = $lw*8;
  } else {
    $ax11 = $rcx+$innerringx;
    $ay11 = $rcy+$y*$ymult;
    $ax12 = $rcx+$innerringx + ($maxtextlen*$tmult*$av1*$tvchartscale);
    $ay12 = $rcy+$y*$ymult;
    $ax21 = $rcx+$outerringx;
    $ay21 = $rcy+$y*$ymult;
    $ax22 = $rcx+$outerringx + ($maxtextlen*$tmult*$av2*$tvchartscale);
    $ay22 = $rcy+$y*$ymult;
    $linewidth = $lw*4;
  }

  $fgcolor = $color1;
  printline($ax11,$ay11,$ax12,$ay12);
  $fgcolor = $color2;
  printline($ax21,$ay21,$ax22,$ay22);

  $fgcolor = $fc;
  $linewidth = $lw;
}

sub dgridline($$) {
  my $v = shift;
  my $chartnum = shift;
  my $fc = $fgcolor;
  my $r;

  if($chartnum == 2) {
    $r = ($maxx*$xmult)+($maxtextlen*$tmult*(1+($tvchartscale*(1+($v/$maxtaxval2)))))+$toffx;
  } else {
    $r = ($maxx*$xmult)+($maxtextlen*$tmult*(1+($tvchartscale*$v/$maxtaxval1)))+$toffx;
  }
  my $ay1 = $rcy;  
  my $ay2 = $rcy+$numtax*$ymult;

  $fgcolor = "#CCCCCC";

  if($plotstyle==2) {
    printcircle($rcx, $rcy, $r, $linewidth);
  } else {
    printline($rcx+$r, $ay1, $rcx+$r, $ay2);
  }

  $fgcolor = $fc;
}

sub dlabel($$) {
  my $y = shift;
  my $str = shift;
  my $clr = shift;
  my $ax;
  my $ay;
  my $fc = $fgcolor;

  $fontsize = $labfontsize;
  
  if(length($clr)>0) {
    $fgcolor = $clr;
  }

  if($plotstyle==2) {
    $ax = radialx($maxx*$xmult+$maxtextlen*$tmult*1.3*$labeloffsetmult, $y);
    $ay = radialy($maxx*$xmult+$maxtextlen*$tmult*1.3*$labeloffsetmult, $y);
    #Print label with the middle aligned closely to the given taxon
    printtext($ax + ((0.90 * sin(deg2rad($y*$ymult))-1) * (length($str)*$tmult*$labelfontmult)/2), $ay, $str);
  } else {
    $ax = $rcx+$maxx*$xmult+$maxtextlen*$tmult*1.3*$labeloffsetmult;
    $ay = $rcy+$y*$ymult;
    printtext($ax, $ay, $str);
  }

  $fgcolor = $fc;

}

sub printline($$$$) {
  my $x1 = shift;
  my $y1 = shift;
  my $x2 = shift;
  my $y2 = shift;
  
  $svglayerstr{$layer} .= "\t<path fill=\"none\" stroke=\"$fgcolor\" stroke-width=\"$linewidth\" d=\"M $x1 $y1 L $x2 $y2\"/>\n";
  
}

sub printarc($$$$$$) {
  my $x1 = shift;
  my $y1 = shift;
  my $x2 = shift;
  my $y2 = shift;
  my $r = shift;
  my $inout = shift;
  
  $svglayerstr{$layer} .= "\t<path fill=\"none\" stroke=\"$fgcolor\" stroke-width=\"$linewidth\" d=\"M $x1 $y1 A $r $r 0 $inout 0 $x2 $y2\"/>\n";
  
}

sub printbox($$$$) {
  my $x1 = shift;
  my $y1 = shift;
  my $x2 = shift;
  my $y2 = shift;
  
  $svglayerstr{$layer} .= "\t<path fill=\"$fgcolor\" stroke=\"none\" stroke-width=\"0\" d=\"M $x1 $y1 L $x1 $y2 L $x2 $y2 L $x2 $y1 Z\"/>\n";
  
}

sub printtext($$$$) {
  my $x = shift;
  my $y = shift;
  my $str = shift;
  my $alignx = shift;
  my $aligny = shift;
  my $alignxstr;
  my $alignystr;

  if($alignx eq "center") {
    $alignxstr = "text-anchor=\"middle\"";
  } else {
    $alignxstr = "";
  }

  if($aligny eq "bottom") {
    $alignystr = "";
  } else {
    $alignystr = "dominant-baseline=\"central\"";
  }

  $svglayerstr{$layer} .= "<text x=\"$x\" y=\"$y\" font-size=\"$fontsize\" font-family=\"$fontfam\" fill=\"$fgcolor\" $alignxstr $alignystr >".$str."</text>\n";
}

sub printtextrot($$$$) {
  my $x = shift;
  my $y = shift;
  my $a = shift;
  my $str = shift;

  $svglayerstr{$layer} .= "<g transform=\"translate($x,$y)\">\n";
  $svglayerstr{$layer} .= "<g transform=\"rotate($a)\">\n";
  $svglayerstr{$layer} .= "<text x=\"0\" y=\"0\" font-size=\"$fontsize\" font-family=\"$fontfam\" fill=\"$fgcolor\" dominant-baseline=\"central\">";
  $svglayerstr{$layer} .= $str;
  $svglayerstr{$layer} .= "</text>\n</g>\n</g>\n";
  
}

sub printcircle($$$$) {
  my $x = shift;
  my $y = shift;
  my $r = shift;
  my $w = shift;

  $svglayerstr{$layer} .= "<circle cx=\"$x\" cy=\"$y\" r=\"$r\" stroke=\"$fgcolor\" fill=\"none\" stroke-width=\"$w\"/>\n";

}

sub printfillcircle($$$) {
  my $x = shift;
  my $y = shift;
  my $r = shift;

  $svglayerstr{$layer} .= "<circle cx=\"$x\" cy=\"$y\" r=\"$r\" stroke=\"$fgcolor\" fill=\"$fgcolor\" stroke-width=\"1\"/>\n";

}

sub printheader() {
  print <<HEAD;
<?xml version="1.0" standalone="no"?>
<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">
<svg width="100%" height="100%" viewBox="0 -1000 $viewwidth $viewheight" preserveAspectRatio="xMinYMin meet" version="1.1" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink">
HEAD
}

sub printfooter() {
  print "</svg>\n";
}

sub findcommonassignment($$) {
  my $liststr1 = shift;
  my $liststr2 = shift;
  my @list1 = split(/;/, $liststr1);
  my @list2 = split(/;/, $liststr2);
  my @commonlist;
  my $i;

  for($i=0; $i<scalar(@list1) && $i<scalar(@list2); $i++) {
    if($list1[$i] eq $list2[$i]) {
      push(@commonlist, $list1[$i]);
    } else {
      last;
    }
  }

  return join(";", @commonlist);
}
