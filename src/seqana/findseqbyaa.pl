#!/usr/bin/perl -w

use seqIterator;
use bioseq;

$infile="/data/genome/hs/hsgngpcrprt";
$seqreader = new seqIterator($infile);

#%posr = (98, G, 167, Y, 168, V, 243, I); 
#%posr = (97, G, 166, Y, 167, V, 242, I); 

while ($seq = $seqreader->next) {
	$header = $seqreader->header;
	#print $seq;
	#$r97 = substr($seq, 97, 1);
	#if ($r97 eq 'G') {
	#	$r166 = substr($seq, 166, 1);
	#	$r167 = substr($seq, 167, 1);
	#	$r242 = substr($seq, 242, 1);
	#	print $header, " 98: ", $r97, " 167: $r166 168: $r167 243: $r242\n";
	#}
	if ($seq =~ /^\w{96,98}([GED])\w{68,70}([YF][VM])\w{74,76}([IV])/) {
		print "$header\n";
		printSeq($seq, *STDOUT);
		print "$1  $2  $3\n";
	}
}
