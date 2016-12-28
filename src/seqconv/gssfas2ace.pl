#!/usr/local/bin/perl -w 
#adopted to convert gss fasta file to ace
#converts fasta sequences to ace format for addition into the database
#only makes a sequence object and a DNA object
#
#usage gssfas2ace fastafile [genomic or cDNA]
#output file will be constructed automatically

$i = 0;
if (@ARGV == 0) {
	print "Usage: fas2ace -t [Type of DNA, genomic or cDNA]";
	print "\n-r Remark -s Segment_of fastaFile\n";
	die;
}
$type = ""; #g, c, or other string
$remark = "";
$segof = "";

while ($ARGV[$i]) {
	if ($ARGV[$i] eq "-t") {
		$i++;
		if ($ARGV[$i] eq "g" || $ARGV[$i] eq "genomic") {
			$type = "Genomic_canonical";
		}
		elsif ($ARGV[$i] eq "c" || $ARGV[$i] eq "cDNA") {
			$type = "cDNA_EST";
		}
		else {$type = $ARGV[$i];}
	}
	elsif ($ARGV[$i] eq "-r") { 
		$remark = $ARGV[++$i]; 
	}
	elsif ($ARGV[$i] eq "-s") { 
		$segof = $ARGV[++$i];
	}
	else { $infile = $ARGV[$i]; }
	$i++;
}

#############
open(INF, $infile) or die "Can't opn $infile: $!\n";
$outf = nameOutf($infile, "ace");
open(OUTF, ">$outf");

%segment=();
$segCnt=0;

$_ = <INF>;
#assumes that the first line is >seqname etc
$seqCount=0;
while ($_ && /^>/) {
	$seqCount++;
	if ($seqCount%100 == 0) {
		print "$seqCount sequences processed\n";
	}

	/clone (\w+)/;
	$seqName = $1;
	$cosmid = substr($seqName, 0, 6);

	$_ = <INF>;	
	$seq = "";
	while ($_ && !/^>/) {
		chomp;
		if (/-/) {
			die "- detected in sequence\n";
		}
		s/\W+//;  #pick sequence only avoiding ^M thing
		$seq .= $_;
		$_ = <INF>;
	}

	$len = length($seq);
	print OUTF "Sequence $seqName\nDNA $seqName  $len\n";
	print OUTF "Segment_of $cosmid\nShotgun\n";

	print OUTF "\nDNA $seqName\n";
	printSeq($seq);
	print OUTF "\n";
	if (!$segment{$cosmid})  {
		$segment{$cosmid}++;
		$segCnt++;
	}
}
foreach $seg (keys %segment) {
	print OUTF "\nSegmentedSeq $seg\nGSS shotgunned\n";
	#assuming all these are from shotgun project
	#you need to check the file and make sure
}
	
print "$segCnt cosmids processed\n";
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

