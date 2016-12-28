### perl ###

use strict;
use Bioseq;

my ($infile, $outfile);
my $i=0;
while ($ARGV[$i]) {
   if ($ARGV[$i] eq "-o") {
      $outfile=$ARGV[++$i];
   }
   else {
      $infile=$ARGV[$i];
   }
   ++$i;
}
if (!$outfile) {
   if ($infile =~ /(.+?)\.fas$/) {
      $outfile = $1 . "_good.fas";
   }
   else {
      $outfile = $infile . ".fas";
   }
}

# first line >seqname
# second line very long not line-breaked
open IN, "<$infile" or die $!;
open OU, ">$outfile" or die $1;

my $title=<IN>;
my $seq=<IN>;
while ($title && $title =~ /^>/) {
   chomp $seq;
   print OU $title;
   printSeq($seq, \*OU);
   $title=<IN>;
   last if (!$title);
   $seq=<IN>;
}

