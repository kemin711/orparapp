### perl ###
##!/usr/bin/perl -w

# merge two fasta file by identifier

use strict;
my $i=0;
my @files=();
while ($ARGV[$i]) {
   push @files, $ARGV[$i];
   ++$i;
}

my %ids=();
my $cnt=0;
foreach my $f (@files) {
   readSeq($f);
}

sub readSeq {
   my $file=shift;
   open IN, "<$file" or die $!;
   $_=<IN>;
   while ($_) {
      while ($_ && !/^>(.+?)\b/) { 
         $_=<IN>; 
      }
      /^>(.+?)\b/;
      my $id=$1;
      ++$cnt;
      if (exists $ids{$id}) {
         $_=<IN>;
         while ($_ && !/^>/) { $_=<IN>; }
      }
      else {
         ++$ids{$id};
         print;
         $_=<IN>;
         while ($_ && !/^>/) { 
            print;
            $_=<IN>; 
         }
      }
   }
}

print STDERR "$cnt sequence, ", scalar(keys %ids), " unique\n";
