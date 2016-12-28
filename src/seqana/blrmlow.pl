### perl ###
##!/usr/bin/perl -w 
#this program remove low-scoring blast hits from a blast results file
#if the first hit is lower than the cutoff score then all matches will
#be discared
#usage blrmlow [inputfile] [cutoff score=optional]
#file: blrmlow
#defined in arguments $cutoff = 40;
#change cutoff value to change score you want to discard if smaller than
if ($ARGV[0]) {
	$outf = $ARGV[0] . ".hi";
}
else {
	print "Usage: blrmlow infile cutoff_score";
	die;
}
if ($ARGV[1]) {
	$cutoff = $ARGV[1];
}
else {
	print "you did not sepecify a cutoff value for Score.  I will use 40\n";
	$cutoff = 40; #default value
}
open(OUT, ">$outf");
$seqCount=0;
$goodseq=0;
$_ = <>;
while (!eof && /^BLAST/) {
	$buff = $_;
	$seqCount++;
	if ($seqCount%200 == 0) {
		print "processing sequence #$seqCount\n";
	}
	$_ = <>;
	while (!/^Sequences producing/ && !/ No hits found / 
	&& !/^Searching\.+done  Database: nr/) {
		$buff .= $_;
		$_ = <>;
	}
	if (/^Sequences producing/) {
	   $ss = index($_, "(bits"); #ss is the score start index
		$buff .= $_;
		$_ =<>;
		$buff .= $_;
		$_ =<>;
		$se = substr($_, $ss);
		$se =~ s/^ +//; #remove leading space
		@sce = split /  /, $se;
		$score = $sce[0];
		#$sce[1] is the E value
		if ($score > $cutoff) {
			print OUT $buff;
			while (!eof && !/^BLAST/) {
				print OUT $_;
				$_ = <>;
			}
			$goodseq++;
		}
		else {
			while (!eof && !/^BLAST/) {
				$_ = <>;
			}
		}
	}
	else {  # those with no hit
		while (!/^BLAST/ && !eof) {
			$_ = <>;
		}
	}
}
print "Total $seqCount sequences ";
print "\n$goodseq have scores higher than $cutoff\n";
