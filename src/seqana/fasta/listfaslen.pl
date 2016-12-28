### perl ###

use FastaReader;
use Fastaseq;
use warnings;
use strict;

#my $fasfile = "Final_Output.fasta";
#my $pat="*fas";
my ($fasfile, $pat, $outfile);

my $i=0;
while ($ARGV[$i]) {
   if ($ARGV[$i] eq '-i') {
      $fasfile = $ARGV[++$i];
   }
   elsif ($ARGV[$i] eq '-p') {
      $pat=$ARGV[++$i];
   }
   elsif ($ARGV[$i] eq '-o') {
      $outfile = $ARGV[++$i];
   }
   else {
      $fasfile = $ARGV[$i];
   }
   ++$i;
}

if (!$fasfile && !$pat) {
   die "You must provide a glob or a file name\n";
}

if ($fasfile) {
   if (!$outfile) {
      $outfile = makeOutputFileName($fasfile);
   }
   open OU, ">$outfile" or die $!;
   processOnefile($fasfile, \*OU);
}
else { #if ($pat) {
   my @infiles = glob($pat);
   if (!$outfile) {
      $outfile = "sequence_length.tab";
   }
   open OU, ">$outfile" or die $!;
   foreach my $f (@infiles) {
      print "processing $f\n";
      processOnefile($f, \*OU);
   }
}

sub makeOutputFileName {
   my $infile = shift;
   my $oufname;
   if ($infile =~ /(.+)\.(fas|fasta)$/) {
      $oufname = $1 . "_length.tab";
   }
   else {
      $oufname = $fasfile . "_length.tab";
   }
   return $oufname;
}

sub processOnefile {
   my $inf=shift;
   my $fho=shift; # output file handle
   my $fasreader = new FastaReader($inf);
   while (my $faseq = $fasreader->next) {
      print $fho $faseq->name, "\t", $faseq->length, "\n";
   }
}
