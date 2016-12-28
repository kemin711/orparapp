#!/usr/local/bin/perl -w

use strict;

# give a prefix
my $prefix = $ARGV[0];
if (!$prefix) {
   $prefix = "bar";
}

my @files=glob("IonXpress*.bam");
my @nomatch = glob("nomatch*.bam");

#print @files, "\n", @nomatch, "\n";

foreach my $f (@files) {
   my $barcode = substr($f, 10, 3);
   print $barcode, "\n";
   my $newfile = $prefix . $barcode . ".bam";
   rename $f, $newfile;
}
if (-e $nomatch[0]) {
   rename $nomatch[0], "nomatch.bam";
}
