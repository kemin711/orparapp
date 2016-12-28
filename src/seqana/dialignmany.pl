### perl ###
##!/usr/bin/perl -w

use Bioseq;

@files = `ls *.pep`;
$end = @files;
$begin=0;

$i=0;
while ($ARGV[$i]) {
	if ($ARGV[$i] eq "-b") { $begin = $ARGV[++$i]; }
	elsif ($ARGV[$i] eq "-e") { $end = $ARGV[++$i]; }
	$i++;
}

# @files = glob("*.pep");   // will not work on some sun's machines
chomp @files;

for ($i=$begin; $i<$end; $i++) {
	$f = $files[$i];
	print STDERR "Aligning $f        # $i .....\n";
	system("dialign -cw $f");
	#print "done return value from & is $rc\n";
}

################################################################
sub cleanDialign {
	my $infile = shift;
	my $oufile = $infile;
	my $id=$infile;
	$infile =~ s/\.pep/\.ali/;
	$oufile =~ s/\.pep/\.al/;
	open IN, "<$infile" or die "Cannot open $infile $!\n";
	open OU, ">$oufile" or die "Cannot open $oufile $~\n";

	print OU "Sameas ID: $id\n";

	$_ = <IN>;
	while (!/Aligned sequences/) { $_ = <IN>; }
	s/^\s+//;
	print OU;
	$_ = <IN>;   # ===========
	$_ = <IN>;   # empty line
	$_ = <IN>;   # listing of sequence names
	while (!/^$/) {  # reading Aligned Sequences:   length:
		print OU; $_ = <IN>;
	}
	print OU "\n";
	while (/^\s/) { $_ = <IN>; }   # upto the actual alignment
	                               # discard all useless lines
	while (!/Alignment \(FASTA format\):/) {
		while (!/^\s/) {    # the alignment lines           
			s/([A-Ya-y])\s([A-Ya-y])/$1$2/g;
			s/- -/--/g;
			print OU;
			$_ = <IN>;
		}
		print OU "\n";
		while (/^\s/ && !/Alignment \(FASTA format\):/) { $_ = <IN>; }
	}

	close IN;
	close OU;
}

=for delete
sub findRepeat {
	my $infile = shift;
	# print "Looking for repeat in $infile ...\n";
	open SIN, "<$infile" or die "cannot open $infile\n";
	my $seq = "";
	$_ = <SIN>;           # >seqname line
	my $seqName = $_;
	while ($_ && /^>/) {
		$_ = <SIN>;          # first line of sequence
		$seq = "";
		while ($_ && !/^>/ ) {
			chomp;
			$seq .= $_;
			$_ = <SIN>;
		}
		
		#@seqchar=('L', 'S', 'T', 'D', 'E', 'K', 'H', 'Q', 'N', 'A', 'G', 
		#           'P', 'Y', 'C', 'W', 'I', 'V', 'M', 'R', 'F');
		#for ($i=0; $i<@seqchar; $i++) {
		#	my $ch = $seqchar[$i];
		#	if ($seq =~ /($ch{8,})/i) {
		#		print "\nCluster: $infile\n$seqName has repeats $1 \n";
		#		DNAseq::printSeq($seq, *STDOUT);
		#		last;
		#	}
		#}
		if ($seq =~ /((.)\2{9,})/ ) {
			print REP "\nCluster: $infile\n$seqName has repeats $1 \n";
			DNAseq::printSeq($seq, *REP);
		}
		if ($seq =~ /((WMDG){3,})/i) {
				print REP "\nCluster: $infile\n$seqName has repeats $1 \n";
				DNAseq::printSeq($seq, *REP);
		}
	}
	close SIN;
}
=cut
