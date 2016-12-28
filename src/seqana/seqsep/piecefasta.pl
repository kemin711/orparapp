### perl ###
##!/usr/bin/perl -w

use Bioseq;
use strict;

my $i=0;
my $nf=600;
my $infile;
my $flist;
while ($ARGV[$i]) {
   if ($ARGV[$i] eq '--number-files'
      || $ARGV[$i] eq '-n') {
      $nf=$ARGV[++$i];
   }
   elsif ($ARGV[$i] eq '-l') {
      $flist=$ARGV[++$i];
   }
   else {
      $infile=$ARGV[$i];
   }
   ++$i;
}

if (!$infile) {
   die "Useage: piecefasta --number-files 600 est.fas\n";
}

my $pieces=breakfasta($infile, $nf);
if ($flist) { # write file list file
   #$flist=fileStem($infile) . '.fls';
   print STDERR "generating file list file: $flist\n";
   open OU, ">$flist" or die $!;
   foreach my $f (@$pieces) {
      print OU "$f\n";
   }
}

