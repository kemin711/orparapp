#!/usr/bin/perl -w
#splits a file containing multiple fasta sequences into multiple files
#each file is named the same name as the sequence name with extension .seq
#this program is mainly for feading into DNA* seqMan sequence alignment

$infile = $ARGV[0];
$outdir = $infile . "dir";
mkdir $outdir, 0755;

$count = 0;
$_ = <>;
while (!eof && /^>/) {
	$count++;
	$seqfile = substr($_, 1);
	chomp $seqfile;
	$seqfile =~ s/ /_/g;
	$seqfile .= ".seq";
	print "$seqfile \n"; #show file name
	open SEQ, ">$outdir/$seqfile";
	$_ = <>;
	while (!/^>/ && !eof) {
		print SEQ $_;
		$_ = <>;
	}
	if (eof) { print SEQ $_; }
	close SEQ;
}

print "total squences: $count\n";

