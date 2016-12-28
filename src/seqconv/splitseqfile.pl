#!/usr/bin/perl -w
# file: splitseqfile

# if the input file for blast is big use this one to do the blast
# because the blastall program dies after certain input of sequences!

use integer;
use bioseq;

my $binSize=300;
my $pieces; # total number of files, more useful sometimes

if (@ARGV < 1) { usage(); }
my $i=0;
while ($ARGV[$i]) {
	if ($ARGV[$i] eq "-i") { $inputfile=$ARGV[++$i]; }
	elsif ($ARGV[$i] eq "-s") { $pieces=$ARGV[++$i]; }
	elsif ($ARGV[$i] eq "-b") { $binSize=$ARGV[++$i]; }
	else { 
      $inputfile=$ARGV[$i]; 
   }
	$i++;
}

if ($pieces) { 
	$seqCnt = `grep '>' $inputfile | wc -l`;
    chomp $seqCnt;
	$binSize = $seqCnt/$pieces + 1; 
}
#if ($binSize > 500) { $binSize = 300; }
splitFastaSeqfile($inputfile, $binSize);

#####################################################
sub usage {
	print STDERR <<ENDU;
splitseqfile -i inputfile -s pieces integer
	-b binsize integer
    or 
    give the input file as the argument.

splitseqfile -s 600  input.fas
 will break input.fas into 600 pieces.

ENDU
	die;
}

__END__

=head1 DESCRIPTION

for breaking large fasta files into smaller pieces


