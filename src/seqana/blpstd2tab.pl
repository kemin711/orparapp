### perl ###
##!/usr/bin/perl -w
# got bug in reading multiple segments of alignment
# needs to fix it before using it

=head1 blpstd2tab Parsing standard blastp output to table

=for purpose
This program is intended to convert the text out put from
standard blastp to the table format.  The standard format
is perfect for human inspection, whereas the table format
is good for the computer to process.  It usually takes a
long time to run a large batch of blast.  This program is
intented for catpturing the essential aspects of the blast
result.  Alternatively you can run the blast with m 8 
option and it will produce the table format.  I am not sure
about the exact columns.  This program will capture the 
length of query and target, score, exp, iden_cnt, sim_cnt,
match_len, and gaps.  Where, iden_cnt is the number of 
identical residues, and sim_cnt is the number of similar
residues.

the columns outputed by this program:
query query_length target target_length score expect match_length 
		identical_residue_count similar_residue_count gapped_residues
separated by tab

expect numbers needs to rounded to match the float data type in 
databases.  My experience is try by error.  e-100 as 0

The blast e-100 is illegal in perl should be converted to 1e-100.
This seems to be enough.
=cut

$ZERO_ROUND = 1e-100;

$_ = <>;
while (!/^Query=/) { $_ = <>; }
while (!eof && /^Query/) {
	/Query=\s*(.+)/;
	$query = $1;
	$_ = <>;
	/(\d+)\s+letters/;
	$query_len = $1;

	# ignore summary portion and remember the no hit possibility
	while (!/^>|\*+ No hits found \*+/) { $_=<>; }
	#***** No hits found ******
	if (/No hits found/) {
		while (!eof && !/^Query=/) { $_=<>; }
		if (eof) { last; }
		next;
	}

	while (/^>/) {  # all targets of one query
		getMatchInfo();
		$_=<>;
		while (!eof && !/^>/ && !/^Query=/) { $_=<>; }
		if (eof) { last; }
	}
}

################### subroutines #######################
sub getMatchInfo {
	$target = substr($_,1);
	chomp $target;
	$_=<>;
	/Length =\s*(\d+)/;
	$target_len = $1;
	while (!/Score =/) { $_ = <>; }
	/Score =\s*(\d+)(?:.+?)Expect =\s*(.+$)/;
	$score=$1;
	$expect=$2;
	$_=<>;
	$gap=0;
	/Identities =\s*(\d+)\/(\d+)(?:.+?)Positives =\s*(\d+)\/\d+/;
	$iden_cnt = $1;
	$match_len = $2;
	$sim_cnt = $3;
	if (/Gaps =\s*(\d+)/) { $gap = $1; }
	$expect =~ s/^e/1e/;
	if ($expect < $ZERO_ROUND) { $expect = 0; }
	print "$query\t$query_len\t$target\t$target_len\t$score\t$expect\t$match_len\t$iden_cnt\t$sim_cnt\t$gap\n";
}
