#!/usr/bin/perl -w 
#converts fasta protein sequence to ace format for addition into the database
#only makes a sequence object and a peptide object
#
#usage fasprt2ace fastafile 
#output file will be constructed automatically

if (@ARGV == 0) {
	print "No file specified for input\n";
	print "Usage: fasprt2ace fastafile\n";
	print "Out put will be put into *.ace file\n";
	die;
}
else {
	$outf = $ARGV[0];
	open(INF, $ARGV[0]) or die"Can't opn $ARGV[0]: $!\n";
}
#print "input file name:  $outf \n";

$outf = nameOutf($ARGV[0], "ace");
open(OUTF, ">$outf");

$_ = <INF>;
#assumes that the first line is >seqname etc
$seqCount=0;
while ($_ && /^>/) {
	$seqCount++;
	if ($seqCount%100 == 0) {
		print "$seqCount sequences processed\n";
	}
	s/>//;
	chomp;
	($seqName, $rest) = split / /, $_, 2;
	$_ = <INF>;	
	$seq = "";
	while ($_ && !/^>/) {
		chomp;
		$seq .= $_;
		$_ = <INF>;
	}

	print OUTF "\nProtein $seqName\n";
	if ($rest) {
		print OUTF "Title  \"$rest\"\n";
	}
	$len = length($seq);
	print OUTF "Peptide $seqName  $len\n";
	print OUTF "\nPeptide $seqName\n";
	printSeq($seq);
	print OUTF "\n";
}

print "$seqCount ace-formated output witten to file: $outf\n";

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

sub printSeq {
	my $seq = $_[0];
   my $i=0;
   my $len = 70;
   my $sub = substr($seq, $i, $len);
   my $seqLen = length($seq);
   print OUTF "$sub\n";
   $i += $len;
   while ($i < $seqLen) {
      $sub = substr($seq, $i, $len);
      print OUTF "$sub\n";
      $i += $len;
   }
}



