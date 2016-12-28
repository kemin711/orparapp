#!/usr/bin/perl -w

use Ace;
use bioseq;

# this program is an extension of the Java ExtractEnd program.
# because ExtractEnd consumes all memory.  I have to write this
# one to finish the job.

# the input format is the tail.ext file
#Genomic NT_008583
#>XM_301267_mod 3'-region   94094 5000
#>XM_291746_mod 3'-region   935435 5000 rc
#>XM_301273_mod 3'-region   1008173 5000
#>XM_291747_mod 3'-region   1011294 5000 rc
#>XM_295810_mod 3'-region   2606429 3674 rc

$hsgn = Ace->connect(-host => "localhost", -port => 3501) or die "cannot open acedb\n";

open IN, "<tail.ext" or die "cannot open tail.fas\n";
open OU, ">tail.fas";

$_ = <IN>;
while ($_) {
	chomp;
	$genomic = substr($_, 8);
	print "Genomic sequence: $genomic\n";

	$dna = convertFas2String($hsgn->fetch("Sequence", $genomic)->asDNA);
	$_ = <IN>;
	while ($_ && /^>/) {
		print OU;
		chomp;  
		@arr = split /\t/;
		$seqname = $arr[0];
		@arr = split / /, $arr[1];
		if (@arr == 1) {
			$subseq = substr($dna, $arr[0]);
		}
		else { # @arr > 1
			$subseq = substr($dna, $arr[0], $arr[1]);
			if (@arr == 3 && $arr[2] eq "rc") {
				$subseq = revcomp($subseq);
			}
		}
		print $seqname, "\n";
		printSeq($subseq, *OU);
		$_ = <IN>;
	}
}
print STDERR "All done!\n";
