### perl ###

use strict;
use Log::Log4perl;
use MiraResult;
use FastaseqQual;

Log::Log4perl->init("/home/zhouke/src/proj/seqconv/log.conf");

my $proj; #= "bar9_trim420cut_7493"; 
my $i = 0;
my $show = 1; # shows the quality of each base
my $internal = 0;
my $cutoff = 0;
my ($margin, $width, $flank);
my $allContigs = 0;
while ($ARGV[$i]) {
   if ($ARGV[$i] eq "-p") {
      $proj = $ARGV[++$i];
   }
   elsif ($ARGV[$i] eq "-i") {
      if ($ARGV[$i + 1] =~ /^\d+$/) {
         $margin = $ARGV[++$i];
      }
      $internal = 1;
   }
   elsif ($ARGV[$i] eq "-l") {
      $cutoff = $ARGV[++$i];
   }
   elsif ($ARGV[$i] eq "-w") {
      $width = $ARGV[++$i];
   }
   elsif ($ARGV[$i] eq "-f") {
      $flank = $ARGV[++$i];
   }
   elsif ($ARGV[$i] eq "-N") {
      $show = 0;
   }
   elsif ($ARGV[$i] eq "--show-all") {
      $allContigs = 1;
   }
   else {
      $proj = $ARGV[$i];
   }
   ++$i;
}

my $miraresult = new MiraResult($proj);
#$miraresult->evaluateResult;
#
my $goodContigs = $miraresult->getGoodContigs;
print "Good contigs: ", join(" | ", @$goodContigs), "\n";

if ($allContigs) {
   foreach my $ctg (@$goodContigs) {
      showOne($ctg, $show);
   }
}
else {
   showOne($goodContigs->[0], $show); 
}


#showInternalLow($fastaq);

###### sub routines #########################
#

sub showOne {
   my $ctgName = shift;
   my $showQuality = shift;

   print "looking at $ctgName\n";

   my $fastaseq = $miraresult->getFastaSequenceById($ctgName);
   my $quality = $miraresult->getQualityById($ctgName);
   my $fastaq = new FastaseqQual($fastaseq, $quality);
   if ($width) {
      $fastaq->width($width);
   }
   if ($flank) {
      $fastaq->flank($flank);
   }

   #print join(" | ", @$quality), "\n";
   if ($showQuality) {
      $fastaq->show;
   }
   if ($internal) {
      $fastaq->showInternalLowQuality($margin);
   }
   elsif ($cutoff) {
      $fastaq->showLowQualityBelow($cutoff);
   }
   else {
      $fastaq->showLowQuality();
   }
}


sub getLogFileName {
   return "./testMiraResult.log";
}

__END__

=head1 NAME

miralowqual - examine the low quality region of mira assembly

=head1 SYNOPSIS

 show the minimum for the whole assembly of bar9_trim420cut_7493
   miralowqual bar9_trim420cut_7493

 show the internal minimum for assembly project bar9_trim420cut_7493:
   miralowqual -i bar9_trim420cut_7493

 show all flanking regions for quality score below 65
   miralowqual -l 65 bar9_trim420cut_7493

=head1 DESCRIPTION

A program to examine the low quality regions of a mira assembly.
There are two main modes of usage: 1) get the minimum quality
2) get qualities for a given cutoff.

This program should be launched from the mira run directory for simplicity.
The argument of the program is the "mira assembly project name", such as
bar9_trim420cut_7493.  It is usually the prefix for the *_assembly directory.

The result of this program can be used to feed the fishing program and 
viewing the assembly at particular regions with gap5.

=head2 Options

   -p project name can also be given as the argument for this program.
   -i end margin to be used for find internal minimum low quality scores.
      if option is given a value it will be used otherwise a default value
      of 12 will be used. The minimum search will ignore the margin 
      bases from each end of the assembly.
   -l cutoff for the quality score. The flanking regions below this
      value will be displayed.
      $cutoff = $ARGV[++$i];
   -w the width to show the quality score of the whole assembly. Default
      is 50 nt.
   -f Size of the flanking region to show for the low quality score.
      default is 18.
   -N tell the program to not show the quality scores for the whole
      sequence.

