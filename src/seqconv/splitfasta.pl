### perl ###
###!/usr/bin/perl -w

use strict;
use warnings;

# file: splitfasta
#splits a file containing multiple fasta sequences into multiple files
#each file is named the same name as the sequence name with extension .fas
#this program is mainly for feading into genscan 
#that does not take multipl sequence stored in the same files

#usage: splitfasta infile

my $FILE_EXTENSION = "fas";
my $DIRSIZE=4000;  # DIRSIZE number of files per directory

if (@ARGV < 1) {
    print STDERR "type perldoc splitfasta to get instruction\n";
    exit;
}

my $i=0;
my $infile;
while ($ARGV[$i]) {
	if ($ARGV[$i] eq "-s") { $DIRSIZE=$ARGV[++$i]; }
	else { $infile = $ARGV[$i]; }
	$i++;
}
print STDERR "$DIRSIZE sequences per directory\n";

# $infile = $ARGV[0];
$infile =~ /(.+)\./;
my $outdirBase = $1;
my $outdir="";

my $count = 0;
open IN, "<$infile" or die $!;

$_ = <IN>;
while (!eof && /^>/) {
	if ($count%$DIRSIZE == 0) {
		$outdir=createDir($count/$DIRSIZE);
		print "opened $outdir  for writing sequences \n"; #show file name
	}
	$count++;
	my $seqfile = substr($_, 1); #assuming one word only >line
	chomp $seqfile;
	my @arr = split / /, $seqfile;
	$seqfile = $arr[0];
	$seqfile .= ".$FILE_EXTENSION";
	# print "$seqfile \n"; #show file name
	open SEQ, ">$outdir/$seqfile";
	print SEQ $_;
	$_ = <IN>;
	while (!/^>/ && !eof) {
		print SEQ $_;
		$_ = <IN>;
	}
	if (eof) { print SEQ $_; }
	close SEQ;
}

print STDERR "total squences: $count\n";

sub createDir {
	my $outdir = $outdirBase . $_[0];
	unless (-d $outdir) {
		$outdir .=  "dir";
		mkdir $outdir, 0755;
	}
	# print STDERR "single sequence files will be written to $outdir\n";
	return $outdir;
}

__END__


=head1 NAME

splitfasta - breaks a fasta file into individual files that 
contains only one sequence per file.


=head1 SYNOPSIS

splitfasta -s 4000 prt.fas


=head1 DESCRIPTION

This program splits the input fasta file into multiple
directories with each directory containing NUM files 
one seq per file.

file: splitfasta

splits a file containing multiple fasta sequences into multiple files
each file is named the same name as the sequence name with extension .fas
this program is mainly for feading into genscan 
that does not take multipl sequence stored in the same files

For programs taking one sequence per file, this is not very productive.
The age of large number of sequences don't like such small programs.

=head1 OPTIONS

    -s <integer>  the number of files (one seq per file) per directory

=head1 AUTHOR

Kemin Zhou 
kemin@att.net

=head1 SEE ALSO

splitseqfile - used to break fasta files into smaller pieces for blast input.
