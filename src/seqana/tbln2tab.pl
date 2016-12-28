#!/usr/bin/perl -w

=head1 tbln2tab Parsing standard (with alignments) tblastn output to table

=for purpose
This program is intended to convert the text out put from
standard tblastn to the table format.  The standard format
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

== note positional information will be captured instead of the 
length of eachparameter done in blpstd2tab

the columns outputed by this program:
query query_length target target_length score expect match_length 
		identical_residue_count positive_residue_count gapped_residues
		Qbegin, Qend, Tbegin, Tend
separated by tab

expect numbers needs to rounded to match the float data type in 
databases.  My experience is try by error.  e-100 as 0

The blast e-100 is illegal in perl should be converted to 1e-100.
This seems to be enough.
=cut

$ZERO_ROUND = 1e-100;
# some applications cannot store values smaller thatn 1e-100

$_ = <>;
while (!/^Query=/) { $_ = <>; }
while (!eof && /^Query=/) {
	/Query=\s*(.+)/;
	$query = $1;
	$_ = <>;
	/(\d+)\s+letters/;
	$query_len = $1;

	# ignore summary portion and remember the no hit possibility
	while (!/^>|\*+ No hits found \*+/) { $_=<>; }
	#***** No hits found ******
	# discard the rest and start a new query
	if (/No hits found/) {
		while (!eof && !/^Query=/) { $_=<>; }
		if (eof) { last; }
		next;
	}

	while ($_ && /^>/) {  # all targets of one query
		getMatchInfo();
		while (!eof && !/^>/ && !/^Query=/) { $_=<>; }
	}
}

################### subroutines #######################
# read one alignment segment, need to repeat this one to 
# exhaust all segments
sub getMatchInfo {
	$target = substr($_,1);
	chomp $target;
	$target =~ s/\s.+//;  # remove trailing title
	$_=<>;
	while (!/\s+Length =/) { $_=<>; }
	/Length =\s*(\d+)/;
	$target_len = $1;
	while (!/Score =/) { $_ = <>; }

	# Segment loop
	while ($_ && /^ Score =/) {
		/Score =\s*(\d+)(?:.+?)Expect.*? =\s*(.+$)/;
		$score=$1;
		$expect=$2;
		$_=<>;
		$gap=0;
		/Identities =\s*(\d+)\/(\d+)(?:.+?)Positives =\s*(\d+)\/\d+/;
		$iden_cnt = $1;
		$match_len = $2;
		$sim_cnt = $3;
		if (/Gaps =\s*(\d+)/) { $gap = $1; }
		$expect =~ s/^e/1e/;  # this seemed corrected
		if ($expect < $ZERO_ROUND) { $expect = 0; }
		# now read the actual alignment
		$_ = <>;
		while (!/^Query:/) { $_=<>; }  # no danger of run away
		# at beginning of alignment
		/.+?(\d+).+?(\d+)/;
		$Qbegin = $1;
		$Qend = $2;
		$_=<>;
		$_=<>; # Target line
		/.+?(\d+).+?(\d+)/;
		$Tbegin = $1;
		$Tend = $2;
		# in case there are more segments
		while (!/^Query:/ && !/^ Score/ && !/^>/ && !/^Query=/) { 
			$_ = <>; 
		}
		while ($_ && /^Query:/) {
			/.+?(\d+).+?(\d+)/;
			$Qend = $2;
			$_=<>; $_=<>;
			/.+?(\d+).+?(\d+)/;
			$Tend = $2;
			while ($_ && !/^Query:/ && !/^ Score/ && !/^>/ && !/^Query=/) { 
				$_ = <>; 
			}
		}
		print $query, "\t", $query_len, "\t",
				$target, "\t", $target_len, "\t",
				$score, "\t", $expect, "\t", $match_len, "\t", $iden_cnt, "\t", 
				$sim_cnt, "\t", $gap,
				"\t", $Qbegin, "\t", $Qend, 
				"\t", $Tbegin, "\t", $Tend, "\n";
	} # segment loop
}
