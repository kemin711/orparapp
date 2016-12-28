### perl ###

use strict;
use Log::Log4perl;
use MiraResult;
use FastaseqQual;

Log::Log4perl->init("/home/zhouke/src/proj/seqconv/log.conf");

my $miraresult = new MiraResult("bar9_trim420cut_7493");
$miraresult->evaluateResult;
my $goodContigs = $miraresult->getGoodContigs;
print "Good contigs: ", join(" | ", @$goodContigs), "\n";

my $ctgName = $goodContigs->[0];
my $fastaseq = $miraresult->getFastaSequenceById($ctgName);
my $quality = $miraresult->getQualityById($ctgName);
#print join(" | ", @$quality), "\n";
testfastaqual();

sub getLogFileName {
   return "./testMiraResult.log";
}

sub testfastaqual {
   my $fastaq = new FastaseqQual("bar9trim420cut_7493", "best assembly same as reference", 
         $fastaseq->sequence, $quality);
   $fastaq->printBoth;

   my ($minimumQuality, $minqIdx) = $fastaq->minQuality;
   print "minimum quality is $minimumQuality at ",
      join(", ", @$minqIdx), " of length of seq ",
      $fastaq->length, "\n";

   ($minimumQuality, $minqIdx) = $fastaq->minQualityInternal(12);
   print "minimum internal quality is $minimumQuality at ", 
      join(", ", @$minqIdx), " of length of seq ",
      $fastaq->length, "\n";
   foreach my $ix (@$minqIdx) {
      print " " x 18, "|\n";
      print $fastaq->subseq($ix - 17, $ix + 19)->sequence, "\n";
   }
}
