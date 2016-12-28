### perl ###
##!/usr/bin/perl -w

#use lib "/home/kzhou/perlmod";
use Bioseq;
use strict;

# file: breakFasfile
# break a large fasta file into smaller pieces
# with -s <num> num number of sequences per file
# if -d <dirsize> then, dirsize file will be put into one directory
# The directory is named from DIR1, DIR2, DIR3,  to whatever

my $numseqPerFile=300;
my $numFiles;

my @arr=split /\//, $0;
my $programName = $arr[$#arr];
if (@ARGV<1) {
   usage();
}
my ($numseqPerFile, $numFilePerDir, $numFiles, $infile);
my $i=0;
while ($ARGV[$i]) {
   if ($ARGV[$i] eq '-s') {
        $numseqPerFile = $ARGV[++$i];
   }
	elsif ($ARGV[$i] eq '-d') {
		$numFilePerDir = $ARGV[++$i];
	}
	elsif ($ARGV[$i] eq '--numfiles') {
		$numFiles = $ARGV[++$i];
	}
   elsif ($ARGV[$i] eq '--help') {
      usage();
   }
   else {
        $infile=$ARGV[$i];
   }
   ++$i;
}

if (!$infile) { usage(); }
# splitFastaSeqfile produces all the files
# in the current directory, need it to 
# put into directoreis
if ($numFiles) {
    $numseqPerFile=computeBinSize($infile, $numFiles);
}
print STDERR "$numseqPerFile sequences per file\n";
# return a reference to an array
my $allFiles = splitFastaSeqfile($infile, $numseqPerFile);
# to be picked up by the shell calling this program
print "Breaking up $infile into ",
	scalar(@$allFiles), " files done!\n";
# now pack into directories

if ($numFilePerDir) {
	my $dirs=packdir();
   print scalar(@$dirs);
}
else {
   print scalar(@$allFiles);
}

####################################

sub usage {
    print STDERR "usage: $programName -s <num seq per file> <input file>\ntype perldoc $programName more information\n";
    exit;
}
   
sub computeBinSize {
    my $file=shift;
    my $numf=shift;
    print STDERR "computing bin size, this may take a while ...\n";
    my $totalseq = `grep -c '>' $file`;
    chomp $totalseq;
    #print "$totalseq sequences in file $file\n";
    if ($totalseq =~ /(\d+) /) {;
      $totalseq=$1;
    }
    #print STDERR "after pattern $totalseq sequences in file $file\n";
    print STDERR int($totalseq/$numf) + 1, " sequence per file\n";
    return int($totalseq/$numf) + 1;
}

sub packdir {
	print STDERR "packing $numFilePerDir file to each directory ...\n";
   my $infileStem=$infile;
	if ($infile =~ /(.+)\.fas$/ || $infile =~/(.+)\.fasta$/) {
		$infileStem = $1;
	}
   my @alldirs=();

	use integer;
	$i=0;
	while ($i < @$allFiles) {
		my $d = $i/$numFilePerDir + 1;
		my $dirname = $infileStem . "DIR" . $d;
		mkdir $dirname;
        push @alldirs, $dirname;
		while ($i<@$allFiles && ($i/$numFilePerDir+1) == $d) {
			system("mv", $allFiles->[$i], $dirname);
			++$i;
		}
	}
	print scalar(@alldirs), " directories created\n";
   return \@alldirs;
}
		

__END__

=head1 NAME 

breakFasfile - break a large fasta file into multiple files
               each with fewer sequences

=head1 SYNOPSIS

breakFasfile -s 100 prt.fas

=head1 DESCRIPTION

This program takes a large fasta file with too many sequences
and break up the large file into smaller files with fewer sequences
in each file.  This is usefull for pre-processing before the 
blast search. All the output files will be stored in the current
directory

=head2 Options

    -s <integer>  specifies the number of sequences per file.
    -d <integer>  specify the number of files per directory
        if not given it will not pack the files into directories.
    --numfiles <integer> To break the whole sequence into specified
        number of files.


=head1 SEE Also

splitfasta: break sequences into one sequence per file, it also
organizes files into directories.

=head1 AUTHOR

Kemin Zhou
Orpara.com


