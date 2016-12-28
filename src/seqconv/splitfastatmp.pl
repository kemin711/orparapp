#!/usr/local/bin/perl -w
#only splits the first two sequence and generate tmp1.fas and tmp2.fas

$FILE_EXTENSION = "fas";

$_ = <>;
for ($i=0; $i <2; $i++) {
	$seqfile = "tmp";
	$seqfile .= "$i.$FILE_EXTENSION";
	print "writing $seqfile to current directory\n"; #show file name
	open SEQ, ">$seqfile";
	print SEQ $_;
	$_ = <>;
	while (!/^>/ && !eof) {
		print SEQ $_;
		$_ = <>;
	}
	if (eof) { print SEQ $_; }
	close SEQ;
}


