package Gff;

=head1 NAME

Gff - read, write, and manipulate GFF files

=head1 SYNOPSIS

 $gff = new Gff( $filehandle );
 foreach $line ($gff->GetLines()) { ... }
 $gff->NewAnnot('Seq1', 'blastx-v2', 'homology',
                45, 123, 0.55, '+', 0, 'HBA_HUMAN');

=head1 DESCRIPTION

A B<Gff> object represents a collection of lines from a GFF file.
Each GFF line is a separate object called a B<GffLine>.

This module also supports an ad-hoc extension to Genie for representing
tiled features in a compact fashion.

To generate a GFF output, usage is as follows:

   my $gff = new Gff;
   my $gffLine = $gff->NewAnnot()
   $gffLine->SetName('SEQ1');
   $gffLine->SetSource('blastx-v2');
   $gffLine->SetFeature('homology');
   $gffLine->SetStart('45');
   $gffLine->SetStop('123');
   $gffLine->SetScore('0.55');
   $gffLine->SetStrand('+');
   $gffLine->SetFrame('0');
   $gffLine->SetExtraFeat('HBA_HUMAN');
   $gffLine->SetTiles( [ s1 s2 s3 ... sN ] );
   $gffLine->SetRegion( [ s1 s2 s3 ... sN ] );

OR just,

   $gff->NewAnnot('Seq1', 'blastx-v2', 'homology',
                  45, 123, 0.55, '+', 0, 'HBA_HUMAN');

OR you can specify the parameters by name:

   $gff->NewAnnot( { Name => 'Seq1',
                     Source => 'blastx-v2',
                     Feature => 'homology',
                     Start => 45,
                     Stop => 123,
                     Score => 0.55,
                     Strand => '+',
                     Frame => 0,
                     ExtraFeat => 'HBA_HUMAN' } );

OR you can copy from another GffLine:

   $gff->NewAnnot( $gffLine );

Remove an entry using Delete(), e.g.

   $gff->Delete( $gffLine );

One or more biosequences may be associated with a GFF file as follows:

   $gff->SetSeq('SEQ1','cggctcggattggcgctggatgata');

To print to a Perl $filehandle, use

   $gff->Print($filehandle);

Sometimes you want to just append data to an existing file.
Use the flags NOSEQ and NOHEAD to control what metadata is included.
To only output the features, not the sequence or header, use

   $gff->Print( $fh, $Gff::NOSEQ|$Gff::NOHEAD ) 

To print to stdout do

   $gff->Print();

or

   $gff->Print(*STDOUT);

A GFF file can be read from a perl FileHandle as follows:

   my $filehandle = new FileHandle "...";
   my $gff = new Gff( $filehandle );

It is also possible to read from *multiple* files.  If the first
argument is a reference to an array, then it is treated as an
array of filehandles.  Each filehandle is read.

   my $gff = new Gff( [$fh1, $fh2] );

You can select a subset of acceptable GFF lines by passing Perl             
regular expressions to match attributes.                                    

For example, the following line only reads in lines with a feature name     
that contains the word 'exon'.                                              

   my $gff = new Gff( $filehandle, 
                      { Feature => '.*exon.*' } );              

See L</RESTRICTIONS> for more information.

The actual biosequence is returned by calling GetSeq(name) as in

   my $seq = $gff->GetSeq('SEQ1');

Iterate through all sequences as follows:

   foreach $line ($gff->GetLines()) { 
     print $line->GetName();
     print $line->GetSource();
     print $line->GetFeature();
     print $line->GetStart();
     print $line->GetStop();
     print $line->GetScore();
     print $line->GetStrand();
     print $line->GetFrame();
     print $line->GetExtraFeat();
     print $line->GetTiles();
   }

to return an array of the unique names of all sequences use

   $gff->GetAllNames()

You can restrict the lines returned from GetLines according to the
name of a sequence; thus iterating over all lines for each sequence is
done as:

   foreach $seqName ($gff->GetAllNames()) {
     foreach $line ($gff->GetLines($seqName)) {
        ...
     }
   }


An extension to the GFF format is implemented to allow describing
a region of tiled GFF features.  There are two different types of
such features, called "Tiles" and "Region" confusingly. 

"Tiles" are defined as a set of consecutive features in which
each tile is independent of the others.  For example:

 SEQ1	EMBL	atg	103	105	0.9	+	0
 SEQ1	EMBL	atg	106	108	0.2	+	0
 SEQ1	EMBL	atg	109	111	0.8	+	0
 SEQ1	EMBL	atg	112	114	0.5	+	0

"Regions" are defined as a set of consecutive features that form
a dependent group.  There may be multiple regions for the same position.
"Regions" are used for database hits.  For example:

 SEQ1	blastx	exon	103	105	0.9	+	0	HIT1
 SEQ1	blastx	exon	106	108	0.2	+	0	HIT1
 SEQ1	blastx	exon	106	108	0.8	+	0	HIT2
 SEQ1	blastx	exon	109	111	0.5	+	0	HIT2

The format for specifying a region uses the '##' comment section so that
standard GFF readers don't complain.  The first example above ("Tiles")
is written as

   SEQ1	EMBL	atg	103	114	2.4	+	0	Reg1
   ##Tiles Reg1
   ##Size 3
   ## 0.9 0.2 0.8 0.5
   ##end-Tiles

The "extra" column of a normal GFF line is used as the name of a region.
For example, above, it's 'Reg1'.  The meaning of the score in a tiled
GFF line is left undefined.

The second example above ("Region") is written as

   SEQ1	blastx	exon	103	108	1.1	+	0	HIT1
   SEQ1	blastx	exon	106	111	1.3	+	0	HIT2
   ##Region HIT1
   ##Size 3
   ## 0.9 0.2
   ##end-Region
   ##Region HIT2
   ##Size 3
   ## 0.8 0.5
   ##end-Region

To create "Tiles" like the first example, above, use SetTiles to pass
a reference to an array of scores.  For example:

   my $gffLine = $gff->NewAnnot('SEQ1', 'EMBL', 'atg', 103, 114, 2.4,
                                '+', 0, 'Reg1');
   $gffLine->SetTiles( [ 0.9, 0.2, 0.8, 0.5 ] );
 
To create a "Region" like the second example, above, use SetRegion in
the same way as SetTiles is used.  For example:

   $gffLine = $gff->NewAnnot('SEQ1', 'blastx', 'exon', 103, 108, 1.1,
                             '+', 0, 'HIT1');
   $gffLine->SetTiles( [ 0.9, 0.2 ] );
   $gffLine = $gff->NewAnnot('SEQ1', 'blastx', 'exon', 106, 111, 1.3,
                             '+', 0, 'HIT2');
   $gffLine->SetTiles( [ 0.8, 0.5 ] );

If a tiled region exists associated with the extra column value,
then GetTiles or GetRegion will return a reference to an array of scores.

   @scores = @{$gffLine->GetTiles( 'Reg1' )};
   @scores = @{$gffLine->GetRegion( 'HIT2' )};

A second extension of the GFF format is to allow a special invalid
score represented by 'Inv'.  This might be interpreted as infinity
in a log-prob context.  For example:

 SEQ1	EMBL	atg	109	111	Inv	+	0

Finally, B<NOTE!!!!>: GFF files index sequences from 1, but any
right-minded, self-respecting programmer numbers from 0.  Therefore,
when reading or writing GFF files, the start and stop values are
adjusted accordingly.  For example, given the GFF line:

 SEQ1	EMBL	atg	103	105	0.9	+	0

After reading in this line from a file, GetStart() returns 102.
Furthermore, calling SetStart(199) will return 199 from GetStart() 
but will print 200.  Usually you can just ignore this fact.

=head1 METHODS

=over 4

=item new( [$fh [, $restrictions ] ] )

Instantiate a new Gff object.  If $fh is specified, then read from the
filehandle.  If $fh is a reference to an array, then read from every
filehandle in the array.  If $restrictions is specified, then exclude
any features that do no match the restriction.  See L</RESTRICTIONS>.

It is assumed that the filehandle(s) are open for reading.

=item Print( $fh [, $flags ] )

Print the Gff object to the specified filehandle.  The optional $flags
modifies what is actually printed.  Currently supported flags are
$Gff::NOHEAD to supress the header (indicating version and date)
and $Gff::NOSEQ to supress the output of the raw sequence.

   $gff->Print( *STDOUT );
   $gff->Print($fh, $Gff::NOHEAD | $Gff::NOSEQ );

It is assumed that the filehandle is open for writing.

=item SetSeq( $label, $dna )

Associates the string of nucleotides $dna with the sequence name
$label.

=item NewAnnot( @fields )

Creates a new Gff entry using the values in @fields.  It is
assumed that @fields is an array containing 8 (or 9 with extra
field) values such that the first is the sequence name, the
second is the source, etc.

NewAnnot returns a GffLine object.

=item NewAnnot( %fields )

Like NewAnnot( @fields ) except the attributes are named in
the hash.  See L</ATTRIBUTE NAMES>.  For example,

   $gff->NewAnnot( { Name => 'Seq1',
                     Source => 'blastx-v2',
                     Feature => 'homology',
                     Start => 45,
                     Stop => 123,
                     Score => 0.55,
                     Strand => '+',
                     Frame => 0,
                     ExtraFeat => 'HBA_HUMAN' } );

NewAnnot returns a GffLine object.

=item NewAnnot( $line )

Allows copying a GffLine object.  This is convenient when
working with two Gff objects and copying lines from one
Gff object to the other.

NewAnnot returns a GffLine object.

=item Delete( $line )

Removes $line from the list of GffLine objects.  $line is
a GffLine object.  For example, to move all lines with source
equal to "SRC" from $gff1 to $gff2, do:

  foreach $line ($gff1->GetLines()) {
     if ($line->GetSource() eq 'SRC') {
        $gff2->NewAnnot($line);
        $gff1->Delete($line);
     }
  }

=item GetSeq( $label )

Returns a string of nucleotide sequences associated with $label.

=item GetLines( [ $label ] )

Returns an array of GffLine objects.  If $label is specified, then
only those lines with the sequence name $label are returned.

=item GetAllNames()

Returns an array of the names of the sequences.  Only the names
for which a known nucleotide sequence exists are returned.

=back

=head2 RESTRICTIONS

The new() method allows filtering of an input Gff object by specifying
restrictions that match a certain criteria.  For example, you may
want to read only Gff lines containing a specific function, source,
score, etc.

The optional second parameter to new() is a hash.  The keys of the
hash are the names of the attributes to restrict.  For example,

   my $gff = new Gff( $filehandle, { Feature => '.*exon.*' } );

There is a special case shorthand for a restriction only on the             
feature attribute:                                                          

   my $gff = new Gff( $filehandle, '.*exon.*' );                             

To retrieve only lines with a sequence name of 'SEQ1', a source name        
starting with 'SRC', on the forward strand between
200 and 2000, with scores less than 1.0 do something like the following.

Start, Stop, and Score are handled specially.  An operator E<lt>,
E<lt>=, E<gt>, E<gt>=, == must
precede a number.

   my $gff = new Gff( $filehandle,                                            
                     { Source => 'SRC.*',                                    
                       Name => 'SEQ1',                                       
                       Strand => '+',                                        
                       Start => '>=200',                                     
                       Stop => '< 2000',                                      
                       Score => '< 1.0' });                                   

=head2 ATTRIBUTE NAMES

The attribute names are

=over 4

=item *Name

=item *Source

=item *Feature

=item *Start

=item *Stop

=item *Score

=item *Strand

=item *Frame

=item *ExtraFeat

=back

=head1 COPYRIGHT

Copyright, David Kulp, 1998-1999.  Unrestricted use is granted.

=cut

use FileHandle;
use strict;

package GffLine;

my $attrList = "Name|Source|Feature|Start|Stop|Score|Strand|Frame|ExtraFeat";

sub new {
    my ($class, @fields) = @_;

    if (ref($fields[0]) =~ /HASH/) {
	my $hash = $fields[0];
	@fields = ();

	my @wrong;
	if (@wrong = grep(!/$attrList/,keys %$hash)) {
	    die "Unrecognized attribute(s): ",join(' ',@wrong);
	}

	my $i=0;
	foreach (split(/\|/,$attrList)) {
	    $fields[$i++] = $hash->{$_};
	}
    }
    elsif (ref($fields[0]) eq 'GffLine') {
	# copy constructor
	@fields = @{$fields[0]};
    }

    my $this = [ @fields ];
    bless $this, $class;

    return $this;
}

sub Print {
    my ($this, $fh) = @_;

    print $fh join("\t", @$this[0..2]);
    print $fh "\t", $$this[3]+1, "\t", $$this[4]+1;
    if ($$this[5] == 0) {
	print $fh "\t0\t";
    } else {
	printf $fh "\t%.6f\t",$$this[5];
    }
    print $fh join("\t", @$this[6..7]);
    if ($$this[8]) {
	print $fh "\t",$$this[8];
    }
    print $fh "\n";

    print $fh $this->PrintTiles($$this[9],"tile") if ($$this[9]);
    print $fh $this->PrintTiles($$this[10],"region") if ($$this[10]);
}

sub SetName    { $_[0][0] = $_[1]; }
sub SetSource  { $_[0][1] = $_[1]; }
sub SetFeature { $_[0][2] = $_[1]; }
sub SetStart   { $_[0][3] = $_[1]; }
sub SetStop    { $_[0][4] = $_[1]; }
sub SetScore   { $_[0][5] = $_[1]; }
sub SetStrand  { $_[0][6] = $_[1]; }
sub SetFrame   { $_[0][7] = $_[1]; }
sub SetExtraFeat { $_[0][8] = $_[1]; }
sub SetTiles  { $_[0][9] = $_[1]; }
sub SetRegion { $_[0][10] = $_[1]; }

# Add the GFF extended Tiles format of the form:
##Tiles {ExtraFeat}
##Size num
##s1 s2 s2 ... sN
##end-Tiles
sub PrintTiles  {
    my ($this, $scores, $tileOrRegion) = @_;

    my $label;
    if ($tileOrRegion =~ /^t/i) {
	$label = "Tiles";
    }
    else {
	$label = "Region";
    }

    # preamble
    my $region = "##$label " . $this->GetExtraFeat() . "\n";
    $region .= "##Size ";
    $region .= ($this->GetStop() - $this->GetStart() + 1) / scalar @$scores;
    $region .= "\n";

    # run length encode
    my @counts;
    my @vals;

    foreach (@$scores) {
	if (/Inv/) {
	    unshift(@vals,'Inv');
	    unshift(@counts,1);
	}
	else {
	    my $sc = sprintf "%.4f", $_;
	    if ($sc eq $vals[0]) {
		$counts[0]++;
	    }
	    else {
		unshift(@vals,$sc);
		unshift(@counts,1);
	    }
	}
    }

    # output scores
    my $counter = 1;
    my $scoreStr;
    while (@vals) {
	my ($val, $count) = (pop(@vals), pop(@counts));
	$scoreStr .= $val;
	$scoreStr .= "($count)" if ($count > 1);
	$scoreStr .= ' ';
	$scoreStr .= "\n##" if ($counter++ % 8 == 0);
    }
    $region .= "##" . $scoreStr . "\n";

    # postamble
    $region .= "##end-$label\n";

    return $region;
}

sub GetName    { return $_[0][0] ;}
sub GetSource  { return $_[0][1] ;}
sub GetFeature { return $_[0][2] ;}
sub GetStart   { return $_[0][3] ;}
sub GetStop    { return $_[0][4] ;}
sub GetScore   { return $_[0][5] ;}
sub GetStrand  { return $_[0][6] ;}
sub GetFrame   { return $_[0][7] ;}
sub GetExtraFeat { return $_[0][8] ;}
sub GetTiles  { return $_[0][9] ;}
sub GetRegion { return $_[0][10] ; }


package Gff;

$Gff::NOSEQ = 1;
$Gff::NOHEAD = 2;

sub new {
    my ($class, $input, $restriction) = @_;
    my $this = {};
    bless $this, $class;

    $this->{feature} = {};	# hash of features -> list of GffLines (broken)
    $this->{list} = [];		# array of GffLines
    $this->{seqs} = {};		# hash of label -> DNA
    $this->{tiles} = {};	# hash of arrays of tiled scores

    if ($restriction) {
	if (ref($restriction) !~ /HASH/) {
	    $restriction = { Feature => $restriction };
	}
	
	if (! keys %$restriction) { goto END_RESTRICT; }

	# check that all keys are acceptable
	my @wrong;
	if (@wrong = grep(!/$attrList/,keys %$restriction)) {
	    die "Unrecognized restriction(s): ",join(' ',@wrong);
	}

	# build an anonymous subroutine to check a GffLine
	my $restrictSub = "sub {  1 ";
	my $key;
	foreach $key ( keys %$restriction ) {
	    my $restrict = $restriction->{$key};
	    if ($key =~ /Stop|Start|Score/) {
		my ($op,$val) = ($restrict =~ /^(.[^\d]?)([e\-\.\d]+)/);
		if ($op !~ /[=><]=?/) {
		    die "Bad comparison operator: $op";
		}
		$restrictSub .= ' && ( $_[0]->Get'.$key."() $op $val ) ";
	    }
	    else {
		my $op;
		if ($restrict =~ /^!/) {
		    $op = '!~';
		    $restrict =~ s/^!//;
		}
		else {
		    $op = '=~';
		}
		$restrictSub .= ' && ( $_[0]->Get'.$key."() $op /$restrict/ ) ";
	    }
	}
	$restrictSub .= "; }";
	$this->{restrict} = eval $restrictSub;
    }

    END_RESTRICT:

    if (defined($input)) {
	if (ref($input) eq 'ARRAY') {
	    foreach (@$input) {
		$this->Read($_);
	    }
	}
	else {
	    $this->Read($input);
	}
    }

    return $this;
}

sub Print {
    my ($this, $fh, $flags) = @_;
    
    $fh = *STDOUT unless defined($fh);

    if (!($flags & $Gff::NOHEAD)) {
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) =
	  localtime(time);
	$mon++;
	# the header info...
	print $fh "##gff-version 1\n";
	print $fh "##source-version Neomorphic-Gff 1.0\n";
	print $fh "##date $year-$mon-$mday\n#\n";
    }	

    if (!($flags & $Gff::NOSEQ)) {
	# the sequences
	my $label;
	foreach $label (keys %{$this->{seqs}}) {
	    my $dna = $this->{seqs}{$label};
	    if (length($dna) > 0) {
		$dna =~ s/(.{50})/$1\n\#\#/g;
		print $fh "##DNA $label\n##";
		print $fh $dna . "\n";
		print $fh "##end-DNA\n#\n";
	    }
	}
    }
    #	print $fh "# Note: scores represent -log(p(x[start..stop])|model)\n\n";
    
    # output GFF lines
    my $line;
    foreach $line (@{$this->{list}}) {
	$line->Print($fh);
    }
}

sub SetSeq {
    my ($this, $label, $dna) = @_;

    $this->{seqs}{$label} = $dna;
}

sub Read {
    my ($this, $input) = @_;

    while (<$input>) {
	# handle special cases:
	if (/^\#\#DNA (\S+)/) {
	    my $seq_name = $1;
	    
	    # read the sequence until we hit 'end-DNA'
	    my $seq = <$input>;
	    while ($seq !~ /\#\#end-DNA/ && !$input->eof()) {
		$seq =~ s/^\#\#//;
		$seq =~ s/\s//g;
		$seq =~ tr/a-z/A-Z/;
		$this->{seqs}{$seq_name} .= $seq;
		$seq = <$input>;
	    }
	}

	if (/^\#\#Tiles\s*(\S+)/) {
	    my $reg_name = $1;
	    my $scores = [];
	    $this->{tiles}{$reg_name} = $scores;

	    # read and stuff
	    my $sizeLine = <$input>;
	    die "Bad tiled region, no size" if $sizeLine !~ /^\#\#Size\s+\d+/;
	    my $regLine = <$input>;
	    while ($regLine !~ /\#\#end-Tiles/ && !$input->eof()) {
		$regLine =~ s/^\#\#//;
		my @lineScores = split(/\s+/,$regLine);
		my $score;
		foreach $score (@lineScores) {
		    my $count;
		    next if $score =~ /^\s*$/;
		    if ($score =~ /\((\d+)\)/) {
			$count = $1;
			$score =~ s/\(.*//;
		    } else {
			$count = 1;
		    }
		    foreach ( 1 .. $count ) {
			push(@$scores,$score);
		    }
		}
		$regLine = <$input>;
	    }
	}

	chop;

	# remove comments
	s/\#.*//;
	next if (/^\s*$/);

	# read in the line
	my @fields = split(/\t+/);

	if (@fields < 8) {
	    warn "Skipping Malformed Line: '$_'.\n";
	    next;
	}

	# add to uniq list of seq names
	$this->{seqs}{$fields[0]} = "" if !defined($this->{seqs}{$fields[0]});

	# subtract one from start, stop
	$fields[3]--; $fields[4]--;

	my $line = new GffLine(@fields);

	if ($this->{restrict}) {
	    push(@{$this->{list}}, $line) if &{$this->{restrict}}($line);
	}
	else {
	    push(@{$this->{list}}, $line);
	}
    }

    # if we've read any tiled regions in, then go back and associate lines
    # with tiled regions.
    if (keys %{$this->{tiles}}) {
	my $line;
	foreach $line (@{$this->{list}}) {
	    if ($this->{tiles}{$line->GetExtraFeat()}) {
		$line->SetTiles($this->{tiles}{$line->GetExtraFeat()});
	    }
	}
    }
}

# remove the line from the list containing the specified ref
# create a temp hash of lines to delete and filter the next time
# GetLines() is called.
sub Delete {
    my ($this, $line) = @_;
    $this->{'delete'}{$line} = 1;
}

sub NewAnnot {
    my ($this, @fields) = @_;

    my $line = new GffLine(@fields);

    push(@{$this->{list}}, $line);

    return $line;
}

sub GetSeq {
    my ($this, $name) = @_;
    
    return $this->{seqs}{$name};
}

# returns all the lines with a name that matches $label, or if no $label
# is given then return all lines.
sub GetLines {
    my ($this,$label) = @_;

    # remove any lines that were registered for deletion
    my $rejects = $this->{'delete'};
    if ($rejects) {
	my $i;
	my $lines = $this->{list};
	for ($i=0; $i <= $#$lines; ) {
	    if ($rejects->{$$lines[$i]}) {
		splice(@$lines,$i,1);
	    }
	    else { $i++; }
	}
	delete($this->{'delete'});
    }

    return @{$this->{list}} unless defined($label);

    my @lines;
    my $line;
    foreach $line (@{$this->{list}}) {
	if (!$label || $line->GetName() eq $label) {
	    push @lines, $line;
	}
    }
    return @lines;
}

sub GetAllNames {
    my ($this) = @_;
    return keys %{$this->{seqs}};
}

1;
