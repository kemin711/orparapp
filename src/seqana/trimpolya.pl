#!/usr/bin/perl -w

use seqIterator;
use bioseq;

$file = $ARGV[0];
if (!$file) {
	die "Usage: trimpolya fasta_file\n";
}

#$file = "/data/genome/hs/gpcr_select.dna";

open OU, ">$file.npa";

$seqReader = seqIterator->new($file);
while ($seq = $seqReader->next) { 
	if ($seq =~ /(a{5,})$/) {
		print $seqReader->header, " POLYA tails\n$1\n";
		$alen = length $1;
		$seq = substr($seq, 0, length($seq)-$alen);
	}
	print OU $seqReader->header, "\n";
	printSeq($seq, *OU);
}
