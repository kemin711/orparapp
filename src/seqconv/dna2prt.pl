### perl ###
##!/usr/bin/perl -w
use DNAseq;

#translates DNA into protein

#usage translate -f DNA sequence file
#                -s start -e end of sequence
#                [+|-][1 2 3]  specifies the frame
# not fully implemented yet!
#
$dna = DNAseq::getSeqFromFile($ARGV[0]);
$pep = DNAseq::translate($dna);
print $pep;

sub usage {
   print "translates DNA into protein\n",
      "usage translate -f DNA sequence file -s start -e end of sequence\n",
      " [+|-][1 2 3]  specifies the frame\n";
   exit(1);
}
