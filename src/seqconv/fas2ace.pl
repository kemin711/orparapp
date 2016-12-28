#!/usr/bin/perl -w 
#converts fasta sequences to ace format for addition into the database
#only makes a sequence object and a DNA object
#
#usage fas2ace fastafile [genomic or cDNA]
#output file will be constructed automatically
# adding one more rarely used option for the fugu genome

use bioseq;

$i = 0;
if (@ARGV == 0) {
	print "Usage: fas2ace -t [g or c; Type of DNA, genomic or cDNA]";
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
	elsif ($ARGV[$i] eq "-tag") { 
		$tag = $ARGV[++$i]; 
	}
	elsif ($ARGV[$i] eq "-s") { 
		$segof = $ARGV[++$i];
	}
	elsif ($ARGV[$i] eq '--fugu') {$FUGU = 1; }
	else { $infile = $ARGV[$i]; }
	$i++;
}

#############
open(INF, $infile) or die "Can't opn $infile: $!\n";
print "Working on $infile ...\n";
$outf = nameOutf($infile, "ace");
open(OUTF, ">$outf");

$_ = <INF>;
#assumes that the first line is >seqname etc
$seqCount=0;
while ($_ && /^>/) {
	$seqCount++;
	if ($seqCount%100 == 0) {
		print STDERR ".";
		if ($seqCount%7000 == 0) { 
			print STDERR "$seqCount sequence processed\n"; 
		}
	}
	chomp;
	($seqName, $rest) = split(/ /, $_, 2);
	$seqName = substr($seqName, 1);  # discard the first >
#for Fugu only
	if ($FUGU) {
		$seqName = "S" . substr($seqName, 9);
	}

=for useless
	if ($rest && $rest =~ /^(.+)/) {
		$rest = "";
	}
=cut

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
	if ($remark) {print OUTF "Remark \"$remark\"\n";}
	if ($segof) {print OUTF "Segment_of \"$segof\"\n";}
	if ($type) {print OUTF "$type\n";}
	if ($rest) {print OUTF "Title \"$rest\"\n";}
	if ($tag) { print OUTF "$tag\n"; }

	print OUTF "\nDNA $seqName\n";
	printSeq($seq, *OUTF);
	print OUTF "\n";
}

print "$seqCount ace-formated output witten to file: $outf\n";

=for removed
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
=cut
