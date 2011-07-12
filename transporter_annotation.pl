#!/usr/bin/env perl 

# By Qinghu Ren
# modified by J. Badger 07/05/2011

use Getopt::Std;
use Switch;
getopts("t:h:c:m:s:n:i:u:d:p:z:v:k:");

# i is the library id/metagenome id
# m is tmhmm search file
# t is transporter blast search trans file
# h is transporter hmm search file
# c is cog search btab file
# s is the pep sequence file
# n is the nraa blastp file
# p is the APIS dataset
# z is the metagenome site id
# u is the cluster transporter search file
# d is the search file directory
# v run verbosely
# k don't exclude on keywords

unless ($opt_d) {
    $initialdir = `pwd`;
    chomp $initialdir;
    $initialdir =~ s/\/$//;
}

if ( $opt_s eq "" ) {
    printf("Usage: $0 [options]\n");
    printf(" -i is the library id/metagenome id\n");
    printf(" -m is tmhmm search file\n");
    printf(" -t is transporter blast search trans file\n");
    printf(" -h is transporter hmm search file\n");
    printf(" -c is cog search btab file\n");
    printf(" -s is the pep sequence file\n");
    printf(" -n is the nraa blastp file\n");
    printf(" -l is the log file\n");
    printf(" -p is the APIS dataset\n");
    printf(" -z is the metagenome site id\n");
    printf(" -u is the cluster transporter search file\n");
    printf(" -d is the search file directory\n");
    printf(" -v run verbosely\n");
    printf(" -k don't exclude on keyword\n");
    die;
}

$opt_m = $initialdir . "\/" . $opt_m;
$opt_t = $initialdir . "\/" . $opt_t;
$opt_h = $initialdir . "\/" . $opt_h;
$opt_c = $initialdir . "\/" . $opt_c;
$opt_s = $initialdir . "\/" . $opt_s;
$opt_n = $initialdir . "\/" . $opt_n;
$opt_l = $initialdir . "\/" . $opt_l;
$opt_u = $initialdir . "\/" . $opt_u;

$fid_list  = "$ENV{'TGG'}/PEP/transporter_faa/fid_tc.list";
$pfam_list = "$ENV{'TGG'}/PEP/transporter_faa/pfam_tigrfam_fid.list";
$cog_list  = "$ENV{'TGG'}/PEP/transporter_faa/cog_fid.list";
( %fid_hash, %pfam_hash, %cog_hash, %transporter_hash, %trans_count ) = ();
@false_COG_array = (
    'COG0574', 'COG0178', 'COG1007', 'COG0484', 'COG1007', 'COG0666',
    'COG1076', 'COG1158', 'COG3291', 'COG1008', 'COG0642', 'COG2208',
    'COG1419', 'COG1157', 'COG2608', 'COG0386', 'COG0443', 'COG0125',
    'COG1609', 'COG0347', 'COG3195', 'COG2362', 'COG0670', 'COG1159',
    'COG4154', 'COG0508', 'COG2524', 'COG3271', 'COG3008', 'COG5001',
    'COG2199', 'COG1180', 'COG1196', 'COG0561', 'COG3641', 'COG3711',
    'COG3933', 'COG3412', 'COG5386', 'COG1511', 'COG2905', 'COG4087',
    'COG2200', 'COG2885', 'COG3920', 'COG2206', 'COG0252', 'COG1902',
    'COG2207', 'COG0550', 'COG0334', 'COG0515', 'COG2114', 'COG1106',
    'COG3064', 'COG0536', 'COG1160', 'COG1213', 'COG1778', 'COG3523',
    'COG1530', 'COG1651', 'COG1538', 'COG2010', 'COG2204', 'COG0810',
    'COG4395', 'COG0576', 'COG2132', 'COG5295', 'COG0723', 'COG3501',
    'COG3147', 'COG2226', 'COG0030', 'COG1195', 'COG1651', 'COG3044',
    'COG0625', 'COG2304', 'COG0741', 'COG5010', 'COG5430', 'COG3012',
    'COG0300', 'COG3767', 'COG0419', 'COG0560', 'COG2226', 'COG2885',
    'COG4934', 'COG4109', 'COG1918', 'COG4485', 'COG0024', 'COG0486',
    'COG2836', 'COG0744', 'COG3910', 'COG0517', 'COG4097', 'COG1192',
    'COG0489', 'COG1629', 'COG0154', 'COG3395', 'COG4886', 'COG1887',
    'COG1294', 'COG0455', 'COG0773', 'COG4253', 'COG1257', 'COG1146',
    'COG1245', 'COG1652', 'COG0543', 'COG1018', 'COG1249', 'COG0845',
    'COG5342', 'COG0103', 'COG0713', 'COG3604', 'COG3848', 'COG1278',
    'COG0469', 'COG1017', 'COG3217', 'COG3420'
);

@false_nraa_array = (
    'GTP',                       'chromosome',
    'oxidoreductase',            'inclusion protein',
    'calcineurin',               'ankyrin',
    'calmodulin',                'polyductin',
    'fibrocystin',               'NLI interacting factor',
    'WD-40 repeat',              'WD-repeat',
    'WD repeat',                 'WD40',
    'Tetratricopeptide',         'SLEI family',
    'transcription termination', 'excinuclease',
    'diguanylate cyclase',       'ComB',
    'CG15040',                   'TPR repeat',
    'TPR-repeat',                'gliadin',
    'extensin',                  'normocyte',
    'mucin',                     'VrrB',
    'interaptin',                'retinoblastoma',
    'fizzy',                     'Lissencephaly',
    'PRL1',                      'coatomer',
    'transducin',                'zinc finger',
    'glutaredoxin',              'EF-hand',
    'clumping factor',           'Elongation factor',
    'GCN20',                     'amidotransferase',
    'LIM factor',                'sperm-associated',
    'ribonucleoprotein',         'coronin',
    'RNA helicase',              'phototropin',
    'kinesin',                   'ubiquitin',
    'loricrin',                  'RNase L inhibitor',
    'STI1',                      'GAF',
    'chemotaxis',                'LacI',
    'kinase',                    'ribosylation',
    'DnaJ',                      'tetratricopeptide repeat',
    'OEP16',                     'lipoxygenase',
    'FRAS1',                     'two-component',
    'LysM',                      'esterase',
    'helicase',                  'peptidase',
    'Caltractin',                'dnaJ',
    'ehydrogenase',              'transglycosylase',
    'GrpE',                      'SMC',
    'ribonuclease',              'xcinuclease',
    'methylase',                 'dehydratase',
    'carbamoyltransferase',      'peroxidase',
    'initiation factor',         'transcriptional regulator',
    'Tim44',                     'cytochrome c assembly',
    'cytochrome c biogenesis',   'CBS',
    'hydrolase',                 'integrase',
    'flagellum',                 'ranscription',
    'HSP',                       'grpE',
    'MinD',                      'synthase',
    'hydrogenase',               'reductase',
    'nitrogenase',               'polymerase',
    'starvation',                'Ankyrin',
    'phosphatase',               'Sigma',
    'racemase',                  'partition',
    'flagellar',                 'ParA',
    'proteinase',                'gyrase',
    'FliI',                      'DNA repair',
    'UvrABC',                    'haperone',
    'AidA',                      'serine/threonine-protein kinase',
    'ribosylation',              'DnaJ',
    'tetratricopeptide repeat',  'OEP16',
    'lipoxygenase',              'FRAS1',
    'Caltractin',                'centrin',
    'transposase'
);

@false_cluster_array = ('CAM_CL_1084');

$OID = $opt_i;
open_it( IN, "$fid_list" );
while (<IN>) {
    chop;
    ( $a, $b, $c, $d ) = split(/\t/);
    $d =~ s/\.//g;
    $fid_hash{$d} = $a;    #FID

}
close IN;

open_it( IN, "$pfam_list" );
while (<IN>) {
    chop;
    ( $a, $b, $c, $d, $e, $f ) = split(/\t|\n/);
    if ( $a !~ /B/ ) {
        $pfam_hash{$d}{"tc"}        = $a;
        $pfam_hash{$d}{"subfamily"} = $c;
        $pfam_hash{$d}{"substrate"} = $f;
    }
}
close IN;

open_it( IN, "$cog_list" );

while (<IN>) {
    chop;
    ( $a, $b, $c, $d, $e, @others ) = split(/\t|\n/);
    $cog_hash{$e}{"tc"}        = $a;
    $cog_hash{$e}{"subfamily"} = $c;
    $cog_hash{$e}{"substrate"} = $d;
}
close IN;

#process transportDB blastp file
open_it( IN, $opt_t );
$temp = "";
while (<IN>) {
    chop;
    ( $a, $b, $c, $d ) = split(/\t|\n/);

    if ( $temp ne $a ) {
        if ( $b ne "1A33" ) {
            $transporter_hash{$a}{"tc"} = $b;
        }
    }
}
close IN;

#process pfam/tigrfam hmm search file
open_it( IN, $opt_h );
$temp = "";
while (<IN>) {
    chop;
    @a      = split(/\t/);
    $evalue = pop(@a);
    if ( $temp ne $a[5] ) {
        if ( $pfam_hash{ $a[0] } && $evalue < 0.0001 ) {
            $transporter_hash{ $a[5] }{"pfam"} = $a[0];    #pfam model name
            $transporter_hash{ $a[5] }{"pfam_evalue"} = $evalue;    #evalue

        }
        $temp = $a[5];
    }
}
close IN;

#process COG search file
open_it( IN, $opt_c );
while (<IN>) {
    chop;
    @fields = split(/\t/);
    if ( $cog_hash{ $fields[3] } ) {
        $transporter_hash{ $fields[0] }{"cog"}        = $fields[3];
        $transporter_hash{ $fields[0] }{"cog_evalue"} = $fields[10];
    }
}
close IN;

#process cluster transporter search file
open_it( IN, $opt_u );
while (<IN>) {
    chop;
    ( $a, $b, $c, $d, $e, $f, $g ) = split(/\t|\n/);

    $transporter_hash{$a}{"cluster"}   = $c;
    $transporter_hash{$a}{"family"}    = $e;
    $transporter_hash{$a}{"subfamily"} = $f;
    $transporter_hash{$a}{"substrate"} = $g;
    $transporter_hash{$a}{"quality"}   = $d;
}
close IN;

$false_count = 0;

# load sequences
open_it( IN, $opt_s );

my ($save_input_separator) = $/;
$/ = "\n>";
while ( $file_array = <IN> ) {
    if ( $file_array =~ /^>(.*)/s ) {
        $file_array = $1;
    }
    if ( $file_array !~ /.*?>$/s ) {
        $file_array .= ">";
    }
    my ( $ORF, $seq ) = ( $file_array =~ /([^\s]+).*?\n(.*)>/s );
    $seq =~ s/\s|\n//g;
    $transporter_hash{$ORF}{"seq"} = $seq if $transporter_hash{$ORF};
}
close IN;

$/    = $save_input_separator;
$file = $opt_n;
$temp = "";

# load nraa hits
open_it( IN, $file );

while (<IN>) {
    chomp;
    @a = split(/\t/);
    if ( $temp ne $a[0] ) {
        if ( $transporter_hash{ $a[0] } ) {
            ( $nraa, $species ) = split( /\[|\]/, $a[14] );
            ($nraa) = split( /\^\|\^/, $nraa );
            $species =~ s/\s+/ /g;
            $species =~ s/;/ /g;
            $evalue                                   = pop(@a);    #evalue
            $transporter_hash{ $a[0] }{"nraa"}        = $nraa;
            $transporter_hash{ $a[0] }{"species"}     = $species;
            $transporter_hash{ $a[0] }{"nraa_evalue"} = $evalue;
        }
        $temp = $a[0];
    }
}
close IN;

# load tmhmm hits
open_it( IN, $opt_m );
while (<IN>) {
    chomp;
    $line = $_;
    if ( $line !~ /\/usr\/local\/bin/ ) {
        if ( $line =~ /\#/ ) {

            if ( $line =~ /^\# ([^\s]+) Length: (\d+)\s*/ ) {
                $a = $1;
                if ( $transporter_hash{$a} ) {
                    $transporter_hash{$a}{"length"} = "len=" . $2;    #length
                }
            }
            elsif (
                $line =~ /^\# ([^\s]+) Number of predicted TMHs:\s+(\d+)\s*/ )
            {
                $a = $1;
                if ( $transporter_hash{$a} ) {
                    $transporter_hash{$a}{"tms"} = "PredHel=" . $2;    #TMS
                }
            }
        }
        else {
            if ( $out !~ /Topology/ ) {
                $out .= "Topology=";
            }

            ( $a, $b, $property, $topology ) = split( /\t/, $line );
            if ( $a ne $temp ) {
                $transporter_hash{$temp}{"topology"} = $out
                  if ( $transporter_hash{$temp} );
                $out  = "Topology=";
                $temp = $a;
            }

            if ( $property eq "inside" ) {
                $out .= "i";
            }
            elsif ( $property eq "outside" ) {
                $out .= "o";
            }
            else {
                $topology =~ /\s+(\d+)\s+(\d+)/;
                $out .= "$1-$2";
            }
        }
    }
}

open_it( IN, $opt_t );
while (<IN>) {
    chop;
    ( $a, $b, $c, $d ) = split(/\t|\n/);
    if ( $transporter_hash{$a} ) {
        $trans_count{$a} = 0 if ( !$trans_count{$a} );
        $trans_count{$a} += 1;
    }
}
close IN;

foreach $transporter ( keys %transporter_hash ) {
    ( $nraa, $species, $evalue, $length, $tms, $topology, $seq, $TC, $cluster )
      = "";

    if ( grep { $_ eq $transporter_hash{$transporter}{"cluster"} }
        @false_cluster_array )
    {
        printf STDERR (
            "Deleting cluster support for %s because it is in bad cluster %s\n",
            $transporter, $transporter_hash{$transporter}{"cluster"}
        );
        delete $transporter_hash{$transporter};
        next;
    }

    $trans_count{$transporter} = 0 if !$trans_count{$transporter};
    if ( grep { $_ eq $transporter_hash{$transporter}{"cog"} }
        @false_COG_array )
    {
        printf STDERR (
            "Deleting cog support for %s because it is in bad cog %s\n",
            $transporter, $transporter_hash{$transporter}{"cog"}
        );
        delete $transporter_hash{$transporter};
        next;
    }

    if ( $cog = $transporter_hash{$transporter}{"cog"} ) {
        $TC        = $cog_hash{$cog}{"tc"};
        $subfamily = $cog_hash{$cog}{"subfamily"};
        $substrate = $cog_hash{$cog}{"substrate"};
    }
    elsif ( $pfam = $transporter_hash{$transporter}{"pfam"} ) {
        $TC        = $pfam_hash{$pfam}{"tc"};
        $subfamily = $pfam_hash{$pfam}{"subfamily"};
        $substrate = $pfam_hash{$pfam}{"substrate"};
        if ( $substrate eq "" && $transporter_hash{$transporter}{"cluster"} ) {
            $substrate = $transporter_hash{$transporter}{"substrate"};
        }
    }
    elsif ( $cluster = $transporter_hash{$transporter}{"cluster"} ) {
        $family    = $transporter_hash{$transporter}{"family"};
        $subfamily = $transporter_hash{$transporter}{"subfamily"};
        $substrate = $transporter_hash{$transporter}{"substrate"};
    }
    else {
        $TC        = $transporter_hash{$transporter}{"tc"};
        $subfamily = "";
        $substrate = "";
    }

    $count2++;

    if ( !$opt_k == 1 ) {

        #check for false hit from nraa top hits
        $false_flag = 0;
        foreach $false_key_word (@false_nraa_array) {
            if ( $transporter_hash{$transporter}{"nraa"} =~ /$false_key_word/ )
            {
                printf STDERR (
"Deleting nraa support for %s because it hits false keyword %s\n",
                    $transporter, $false_key_word
                );
                $false_flag = 1;
            }
        }
    }

    if ( $false_flag == 1 ) {
        $false_count++;
        next;
    }

    if ( $substrate eq "" ) {
        $substrate = get_substrate_by_TC($TC);
    }

    $count_evidence = 0;
    if ( $trans_count{$transporter} > 0 ) { $count_evidence++ }
    if ( $transporter_hash{$transporter}{"pfam"} ne "" ) { $count_evidence++ }
    if ( $transporter_hash{$transporter}{"cog"}  ne "" ) { $count_evidence++ }
    if ( $transporter_hash{$transporter}{"cluster"} ne "" ) {
        $count_evidence++;
    }

    if ( $family eq "" ) {
        $family = $fid_hash{$TC};
    }

    if ( $family =~ /\// && $TC ne "" ) {
        $family = $fid_hash{$TC};
    }

    if ( $subfamily eq "Enzyme I" || $subfamily eq "Hpr" ) {
        $family = "GPTS";
    }

    print $transporter 
      . $OID . "\t" 
      . $opt_p . "\t" 
      . $opt_z . "\t" 
      . $family . "\t"
      . $subfamily . "\t"
      . $count2 . "\t"
      . $substrate . "\t"
      . $transporter_hash{$transporter}{"tc"} . "\t"
      . $trans_count{$transporter} . "\t"
      . $transporter_hash{$transporter}{"pfam"} . "\t"
      . $transporter_hash{$transporter}{"pfam_evalue"} . "\t"
      . $transporter_hash{$transporter}{"cog"} . "\t"
      . $transporter_hash{$transporter}{"cog_evalue"} . "\t"
      . $transporter_hash{$transporter}{"cluster"} . "\t"
      . $transporter_hash{$transporter}{"quality"} . "\t"
      . $count_evidence . "\t"
      . $transporter_hash{$transporter}{"length"} . "\t"
      . $transporter_hash{$transporter}{"tms"} . "\t"
      . $transporter_hash{$transporter}{"topology"} . "\t"
      . $transporter_hash{$transporter}{"nraa"} . "\t"
      . $transporter_hash{$transporter}{"species"} . "\t"
      . $transporter_hash{$transporter}{"nraa_evalue"} . "\t"
      . $transporter_hash{$transporter}{"seq"} . "\n";

    $family = "";
}


###########################
#Subfunctions
###########################

# opens file for reading whether it is uncompressed, gzipped, or bzip2ed
sub open_it {
    my ( $handle, $file ) = @_;
    if ( $file =~ /\.gz/ ) {
        open $handle, "gunzip -c $file|", or die "Can't open $file\n";
    }
    elsif ( $file =~ /\.bz2/ ) {
        open $handle, "bzcat $file|", or die "Can't open $file\n";
    }
    else {
        open $handle, $file, or die "Can't open $file\n";
    }
    print STDERR "Opening $file for reading...\n" if ( $opt_v == 1 );
}

sub get_substrate_by_TC($TC) {
    switch ($TC) {
        case "2A12" { $substrate = "ADP:ATP antiporter" }
        case "2A18" { $substrate = "amino acid" }
        case "2A59" { $substrate = "arsenite" }
        case "2A31" { $substrate = "chloride/bicarbonate anion exchange" }
        case "2A25" { $substrate = "sodium ion:alanine symporter" }
        case "2A49" { $substrate = "chloride ion channel" }
        case "2A50" { $substrate = "glycerol:proton symporter " }
        case "2A3"  { $substrate = "amino acid" }
        case "3A10" { $substrate = "proton (vacuolar)" }
        case "3A4"  { $substrate = "arsenite (ArsA)" }
        case "2A45" { $substrate = "arsenite (ArsB)" }
        case "2A28" { $substrate = "sodium ion:bile acid symporter" }
        case "2A46" { $substrate = "benzoate" }
        case "2A19" { $substrate = "proton:calcium ion antiporter" }
        case "2A77" { $substrate = "cadmium ion" }
        case "2A24" { $substrate = "sodium ion:citrate symporter" }
        case "2A4"  { $substrate = "cation efflux" }
        case "2A11" { $substrate = "proton:citrate symporter" }
        case "1A11" { $substrate = "ammonium" }
        case "2A41" { $substrate = "sodium ion:nucleoside symporter" }
        case "2A36" { $substrate = "sodium ion:proton antiporter" }
        case "2A37" { $substrate = "potassium/sodium ion:proton antiporter" }
        case "2A63" { $substrate = "sodium ion:proton antiporter" }
        case "2A61" { $substrate = "C4-dicarboxylate" }
        case "2A13" { $substrate = "C4-dicarboxylate" }
        case "9A11" { $substrate = "copper ion" }
        case "9A12" { $substrate = "copper ion" }
        case "2A7"  { $substrate = "drug/metabolite?" }
        case "1A6"  { $substrate = "sodium ion channel" }
        case "2A57" { $substrate = "nucleoside" }
        case "2A27" { $substrate = "sodium ion:glutamate symporter" }
        case "3A2"  { $substrate = "protons" }
        case "2A71" { $substrate = "folate/biopterin" }
        case "9A8"  { $substrate = "ferrous ion" }
        case "2A44" { $substrate = "formate/nitrite" }
        case "1A10" { $substrate = "glutamate gated ion channel" }
        case "2A8"  { $substrate = "gluconate" }
        case "2A42" { $substrate = "aromatic amino acid" }
        case "2A72" { $substrate = "potassium ion uptake" }
        case "2A14" { $substrate = "L-lactate" }
        case "2A78" { $substrate = "branched-chain amino acid efflux (AzlC)" }
        case "2A26" { $substrate = "branched-chain amino acid" }
        case "2A75" { $substrate = "lysine efflux" }
        case "9A2"  { $substrate = "mercuric ion" }
        case "9A19" { $substrate = "magnesium ion" }
        case "1A8"  { $substrate = "glycerol uptake" }
        case "1A35" { $substrate = "magnesium/cobalt ion" }
        case "1A22" {
            $substrate = "large-conductance mechanosensitive ion channel"
        }
        case "1A23" {
            $substrate = "small-conductance mechanosensitive ion channel"
        }
        case "2A70" { $substrate = "sodium ion:malonate symporter" }
        case "2A54" { $substrate = "tricarboxylates" }
        case "2A39" {
            $substrate = "cytosine/purines/uracil/thiamine/allantoin"
        }
        case "2A40" { $substrate = "xanthine/uracil" }
        case "2A33" { $substrate = "sodium ion:proton antiporter" }
        case "2A34" { $substrate = "sodium ion:proton antiporter" }
        case "2A35" { $substrate = "sodium ion:proton antiporter" }
        case "2A62" { $substrate = "sodium ion:proton antiporter" }
        case "2A52" { $substrate = "nickel ion" }
        case "2A55" { $substrate = "manganese/iron ion" }
        case "2A22" { $substrate = "sodium ion:solute symporter" }
        case "2A60" { $substrate = "organic anion" }
        case "9A10" { $substrate = "iron ion" }
        case "2A67" { $substrate = "oligopeptide" }
        case "2A9"  { $substrate = "60 KD inner membrane protein OxaA homolog" }
        case "9A17" { $substrate = "lead uptake" }
        case "2A20" { $substrate = "phosphate" }
        case "2A58" { $substrate = "sodium ion:phosphate symporter" }
        case "9A4"  { $substrate = "potassium ion uptake?" }
        case "9A40" { $substrate = "heavy metal ion" }
        case "2A17" { $substrate = "proton:dipeptide/tripeptide symporter" }
        case "9A18" { $substrate = "SbmA homolog" }
        case "2A76" { $substrate = "amino acid efflux" }
        case "2A6"  { $substrate = "multidrug efflux" }
        case "1A3"  { $substrate = "calcium ion channel" }
        case "2A21" { $substrate = "sodium ion:proline symporter" }
        case "2A53" { $substrate = "sulfate" }
        case "2A64" { $substrate = "protein export" }
        case "2A56" { $substrate = "C4-dicarboxylate" }
        case "2A38" { $substrate = "potassium ion uptake" }
        case "1A28" { $substrate = "urea" }
        case "1A1"  { $substrate = "potassium ion channel" }
        case "2A5"  { $substrate = "zinc ion" }
        case "2A23" { $substrate = "proton/sodium ion:glutamate symporter" }
        case "2A2"  { $substrate = "proton/sodium ion:sugar symporter" }
        case "2A1"  { $substrate = "multidrug efflux" }
        case "2A51" { $substrate = "chromate ion" }
        case "2A16" { $substrate = "tellurite" }
        case "2A68" { $substrate = "aminobenzoyl-glutamate" }
        case "2A66" { $substrate = "multidrug efflux" }
        case "2A10" { $substrate = "2-keto-3-deoxygluconate" }
        case "2A47" { $substrate = "sodium ion:anion symporter" }
        case "2A15" { $substrate = "glycine betaine/carnitine/choline" }
        case "1A29" { $substrate = "amide/urea" }
        case "2A81" { $substrate = "aspartate:alanine antiporter" }
        case "2A80" { $substrate = "tricarboxylate" }
        case "1A33" { $substrate = "heat shock protein 70 homolog" }
        case "2A83" { $substrate = "sodium ion:bicarbonate symporter" }
        case "2A85" { $substrate = "fusaric acid efflux?" }
        case "2A86" { $substrate = "autoinducer-2 (AI-2) export?" }
        case "3A12" { $substrate = "cell division/septum DNA translocation" }
        case "1A56" { $substrate = "copper ion uptake" }
        case "2A30" { $substrate = "potassium ion:chloride ion symporter" }
        case "1A15" { $substrate = "protein translocase Sec62 homolog" }
        case "1A18" {
            $substrate =
              "protein import-related chloroplast anion-selective channel"
        }
        case "1A2" {
            $substrate = "ATP-sensitive inward rectifier potassium ion channel"
        }
        case "1A47" { $substrate = "nucleotide-sensitive anion channel" }
        case "9A15" { $substrate = "tryptophan uptake" }
        case "1A4" {
            $substrate = "transient receptor potential cation channel "
        }
    }

    return $substrate;
}
