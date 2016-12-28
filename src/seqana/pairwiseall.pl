#!/usr/bin/perl -w

# pairwiseall
# driver to align all sequences in this directory and report the
# score or sequence identities, uses the external align program align

$cpu = 2;
$alignProgram = "align";
$allAlignedFile = "aligned.all";
open OU, ">$allAlignedFile";
open SIM, ">similarity.all";

$file_ext = $ARGV[0];

$i=0;
while ($ARGV[$i]) {
	if ($ARGV[$i] eq '-p') { $alignProgram = $ARGV[++$i]; }
	elsif ($ARGV[$i] eq '-e') { $file_ext = $ARGV[++$i]; }
	elsif ($ARGV[$i] eq '-o') { $option = $ARGV[++$i]; }
	else { $file_ext = $ARGV[$i]; }
	$i++;
}

if (!$file_ext) { $file_ext="fas"; }

@files = glob("*.$file_ext");
if (@files == 0) {
	@files = glob("*.dna");
	if (@files == 0) { @files = glob("*.pep"); }
}
if (@files == 0) {
	die " no file found with extension $file_ext, *.pep, *.dna, or *.fas\n";
}
for ($i=0; $i<@files; $i++) {
	for ($j=$i+1; $j<@files; $j++) {
		print $files[$i], " to ", $files[$j], "\n";
		$cmd = $alignProgram . " " . $files[$i] . " " . $files[$j] . " " . $option;
		$output = `$cmd`;
		getSim($output);
		print OU $output, "\n";
	}
}

# different programs produce different formats
# this only workw with align!
sub getSim {
	my $rslt = shift;
#
#ALIGN calculates a global alignment of two sequences
# version 2.0u Please cite: Myers and Miller, CABIOS (1989) 4:11-17
#AP006557                                           29727 nt vs.
#AP006558                                           29725 nt
#scoring matrix: DNA, gap penalties: -16/-4
#100.0% identity;     Global alignment score: 148569
#

	$rslt =~ /^.+$/mg;  # read one line
	$rslt =~ /^.+$/mg;
	$rslt =~ /^(\w+)\s+(\d+).+$/mg;
	$seq1 = $1;
	$seq1len = $2;
	$rslt =~ /^(\w+)\s+(\d+).+$/mg;
	$seq2 = $1;
	$seq2len = $2;
	$rslt =~ /^.+$/mg;
	$rslt =~ /^(\d+\.\d+).+$/mg;
	$percent = $1;
	print SIM "$seq1\t$seq1len\t$seq2\t$seq2len\t$percent\n";
}
