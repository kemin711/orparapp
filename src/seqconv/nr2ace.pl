#!/usr/bin/perl -w 
#converts nr protein sequences to ace format for addition into the database
#only makes a protein and peptide objects
#
#usage nr2ace fastafile 
#output file will be constructed automatically

if (@ARGV == 0) {
	print "No file specified for input\n";
	print "Usage: fas2ace fastafile \n";
	die;
}
else {
	$outf = $ARGV[0];
	open(INF, $ARGV[0]) or warn "Can't opn $ARGV[0]: $!\n";
}

$outf = nameOutf($ARGV[0], "ace");
open(OUTF, ">$outf");

$ln = <INF>;
#assumes that the first line is >seqname etc
$seqCount=0;
while ($ln && $ln =~ /^>/) {
	$seqCount++;
	#if ($seqCount%1000 == 0) {
	#	print "$seqCount sequences processed\n";
	#}
	$ln =~ s/>//;
	chomp $ln;
	($seqName, $rest) = split / /, $ln, 2;
	if (length($rest) > 200) {
		$rest = substr($rest, 0, 200);
	}
	$rest =~ s/\"/ /g;
	@tmp = split /\|+/, $seqName;
	$seqkey = $tmp[1];
	$seqkey = "nr:" . $seqkey;

	$len=0;
	$sequence = "";
	$ln = <INF>;	
	while ($ln && $ln !~ /^>/) {
		$sequence .= $ln;
		$len += (length($ln)-1);
		$ln = <INF>;
	}
	##########output########333
	print OUTF "\nProtein \"$seqkey\"\n";
	print OUTF "Title  \"$rest\"\n";
	if ($sequence =~ s/U/X/g) {
		print OUTF "DB_remark \"Use X as U\"\n";
	}
	print OUTF "Peptide  \"$seqkey\"  $len\n";
	print OUTF "\nPeptide  \"$seqkey\"\n";
	print OUTF $sequence;
	print OUTF "\n";
}

print "$seqCount ace-formated output witten to file: $outf\n";

##############subroutines################

sub nameOutf {
#(infile, extension), returns outfile with .extension
my ($of);
	$i=-1;
	if (($i=index($_[0], ".")) > -1) {
		$of = substr($_[0], 0, $i);
	}
	else {
		$of = $_[0];
	}
	$of .= ".";
	$of .= $_[1];
	return $of;
}

	

