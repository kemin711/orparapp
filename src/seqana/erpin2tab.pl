#!/usr/bin/perl -w

# a simple program to covnert erpin output into table format 
# for input into a relational database


# typical data input
# Training set:  "polya.epn":
#       2327 sequences of length 206
# Database:   "tail.fas"
# Cutoff:     10.05
# 
# >XM_301267_mod 3'-region   94094 5000
# FW   1    3115..3166     10.67
# AATTTT.CTTTTACCCAGAAATTCTACCTTTGTGATTTTTTTTTTTTTTGAGA
# FW   2    4242..4293     13.48
# AATAAA.TAAAAAGAAAAATAAATAATTAAAAGAAATTAAGTTGGTTCCATCT
# >XM_291746_mod 3'-region   935435 5000 rc
# FW   1    1752..1803     14.16
# AATAAA.ATCCTTTTAAATTTTTTTCCATGTGTAATTTTAGGATGATATATAT
# FW   2    2076..2127     10.46
# AATAAA.AATAGTTCTCCCTCTTTAGTCCATCCAATCCTAAAGAATTTTTATC
# FW   3    2447..2498     13.50
# AATAAA.CATCATATTGTTGGATTGTAAAATTTTTGGATTCTCAAGTACTCAG
# FW   4    4786..4837     10.13
# AAGAAA.AAAAGTCTGTACATGTTTGGTACAGATGCAGTTTTTTCTTAAAAAT

$_ = <>;
while (!/^>/) { $_ = <>; }
while ($_ && $_ ne "\n" && /^>/) {
	@arr = split /\s+/, substr($_, 1);
	$seqname = $arr[0];
	$_ = <>;
	while ($_ && $_ ne "\n" && !/^>/) {
		chomp;
		@arr = split /\s+|\.\./;
		$_ = <>;
		chomp;
		push(@arr, $_);
		print $seqname, "\t", join("\t", @arr), "\n";
		$_ = <>;
	}
}
print STDERR "all done\n";
