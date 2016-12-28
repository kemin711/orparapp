### perl ###
#!/usr/bin/perl -w

use Bioseq;
use Pod::Usage;
use warnings;
use strict;

my $action = "d"; #default delete DNA
# action = one of, exlusive, cannot do combiations
# d: delete, sub: subsequence, r reverse-complement

my $start = 0;  # default: start of sequence 0-based index
my $end = -1;   # default: end of the complete sequence
if (@ARGV == 0) {
   pod2usage();
   exit(1);
   #die "usage seqed -s start -e end -r -a int -sub infile\n";
}
my ($insertString, $infile, $oufile);
my $i=0;
while ($ARGV[$i]) {
	if ($ARGV[$i] eq "-s") { $start = $ARGV[++$i]-1; }
	elsif ($ARGV[$i] eq "-e") {
		$end = $ARGV[$i+1]-1;
		$i++;
	}
	elsif ($ARGV[$i] eq "-r") {$action = "rc";}  #reverse complement
	elsif ($ARGV[$i] eq "-sub") { $action = "sub"; }
	elsif ($ARGV[$i] eq "-a") {
		$action = "add";
		$i++;
		$start = $ARGV[$i++];  #insertion position
		$insertString = $ARGV[$i];
	}
   elsif ($ARGV[$i] eq "-o") { $oufile = $ARGV[++$i]; }
	else {$infile = $ARGV[$i];}
	$i++;
}

my $header;
my $seq = getSeqFromFile($infile, $header);
$seq = uc($seq);  # upper case translation
my $seqlen = length($seq);
print STDERR "length of original sequence: $seqlen\n";
###############################

if ($action eq "rc") {
	$seq = revcomp($seq);
}
elsif ($action eq "add") {
	my $tmp = substr($seq, 0, $start);
	$tmp .= $insertString;
	$tmp .= substr($seq, $start);
	$seq = $tmp;
}
elsif ($action eq "sub") {
	print STDERR "Start $start End $end extracted\n";
	if ($end == -1) { 
		$seq = substr($seq, $start);
	}
	else {
		$seq = substr($seq, $start, $end-$start+1);
	}
}
else {
	if ($start == 0) {
		$seq = substr($seq, $end+1);
	}
	else {
		if ($end == -1) {
			$seq = substr($seq, 0, $start);
		}
		else {
			my $tmp = substr($seq, 0, $start);
			$seq = substr($seq, $end+1);
			$seq = $tmp . $seq;
		}
	}
}
# output
my @arr = split /\s+/, $header;
## length after the operation, $subseqlen
my $subseqlen = length($seq);
if ($action eq "rc") {$arr[0] = $arr[0] . "rc";}

my $title = $arr[0];
$start++;  # converts back to 1-based index
if ($end == -1) { $end = $seqlen; }
else { $end++; }

if ($action eq "sub") { $title .= " subsequence $start $end"; }

if ($oufile) {
   open OU, ">$oufile" or die $!;
   print OU "$title  $subseqlen nt\n";
   printSeq($seq, *OU);
}
else {
   print "$title  $subseqlen nt\n";
   printSeq($seq, *STDOUT);
}

__END__

=head1 NAME

seqed - edit fasta sequences and write result to a new file

=head1 SYNOPSIS

seqed -s start -e end -r -a int -sub -o outfile infile
   -r for reverse complement.
   -sub for taking subsequence.

=head1 DESCRIPTION

This is a simple command line fasta sequence editor.

delets DNA or Protein (not implemented yet) sequences 
index starts at 1

=head2 Options

   -s start of the range 1-based index. default start=1, 
   -e end of the raqnge 1-based index inclusive
   -r to get reverseComplement of this sequence
   -a number string add a string at/after the specified position
   -d -sub -r are actions. d is for deletion, r for reverse complement, 
      and sub for subsequence.
   -o outputfile the name of the output file for the result
