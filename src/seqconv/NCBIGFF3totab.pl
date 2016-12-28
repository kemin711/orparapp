#!/jgi/tools/bin/perl -w

use GFF3Converter;

#my $inputfile="/home/kzhou/work/analysis/externalGenomes/Ustilago_maydis/NC_008368.gff";
#my $inputfile="/home/kzhou/work/analysis/externalGenomes/Ustilago_maydis/NW_100967.gff";
my @infiles=();
my $i=0;
while ($ARGV[$i]) {
    if ($ARGV[$i] eq '-i') { $inputfile=$ARGV[++$i]; }
    else {
        $inputfile=$ARGV[$i];
    }
    ++$i;
}
if ($inputfile) {
    push @infiles, $inputfile;
}
else {
    @infiles = glob("*.gff");
}

my $gf=GFF3Converter->new;
if (@infiles > 0) {
    foreach $f (@infiles) {
        $gf->processFile($f);
    }
    $gf->showOrganismInformation;
    print STDERR "GFF Done\n";
}

$gf->sequenceToTable("fna", "faa");

print STDERR "All done!\n";
