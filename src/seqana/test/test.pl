#!/usr/bin/perl -w

use AlnVariantHelper;

my $tmpfile="bar10_cut180.fastq";
print "$tmpfile has stem: ", getFileStem($tmpfile), "\n";


$rslt =<<ENDS;
ALIGN calculates a global alignment of two sequences
 version 2.0u Please cite: Myers and Miller, CABIOS (1989) 4:11-17
AP006557                                           29727 nt vs.
AP006558                                           29725 nt
scoring matrix: DNA, gap penalties: -16/-4
100.0% identity;     Global alignment score: 148569
ENDS

	$rslt =~ /^.+$/mg;  # read one line
	$rslt =~ /^.+$/mg;
	$rslt =~ /^(\w+)\s+(\d+).+$/mg;
	print "seq1 $1 len: ", $2;
	$rslt =~ /^(\w+)\s+(\d+).+$/mg;
	print " seq2 $1 len:", $2;
	$rslt =~ /^.+$/mg;
	$rslt =~ /^(\d+\.\d+).+$/mg;
	print " percent $1\n";
