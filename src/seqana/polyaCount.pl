### perl ###
##!/usr/bin/perl -w

#file: polyaCount.pl
# count the poly A length of mRNA from all mRNA sequences

use Ace;
use bioseq;

$hsgn = Ace->connect(-host => "localhost", -port => 3501) or die "Cannot open hsgenome on poart 3501\n";
@seq = $hsgn->keyset('mRNA');

open OU, ">polyN.fas";

for ($i=0; $i<@seq; $i++) {
	$dna = convertFas2String($seq[$i]->asDNA);  # in fasta format
	print $seq[$i], "\t";
	if ($dna =~ /(a+)$/) {
		print length($1);
	}
	else { print 0; }
	print "\t";
	# only the first polyN is counted, indicatding more polyN in the sequence
	if ($dna =~ /(n{5,})/) { 
		print length($1); 
		print OU '>', $seq[$i], "\n";
		printSeq($dna, *OU);
	}
	else { print 0; }
	print "\n";
}
