### perl ###
##!/usr/bin/perl -w

use strict;
use warnings;
use FastaReader;
use Fastaseq;

# this is only accessible from Kraken
# to access from europa add prefile: /net/kraken
my $guideFile = "/ng14/zhouke/virology/hcvseq/fdaref/ddlref.note";
my $workdir="/ng14/zhouke/virology/hcvseq/fdaref";
my @subnames=("NS34A", "NS5A", "NS5B", "NS5AB");

if ($ENV{HOSTNAME} eq 'europa') {
   $guideFile = '/net/kraken' . $guideFile;
   $workdir = '/net/kraken' . $workdir;
}

generateSub();

sub readGuide {
   my $gf = shift;
   my %result;
   open IN, "<$gf" or die "$! $gf\n";
   <IN>;
   while (/^#/) { # read comment and discard
      <IN>;
   }
   <IN>;
   # 0 hcvgenotype, 4 NS34a, 5 NS5A, 6 NS5B
   while (<IN>) {
      chomp;
      my @row=split /\t/;
      $result{$row[0]} = [$row[4], $row[5], $row[6] ];
   }
   return \%result;
}

sub writeSub {
   my $seq=shift;
   my $range=shift;
   my $sname=shift; # subnames element
   print "using subname: $sname\n";

   my ($start, $end)=split /-/, $range;
   my $sub=$seq->subseq($start, $end);
   #print "old subseqname: ", $sub->name, "\n";
   $sub->name($seq->name . $sname);
   #print "new subseqname: ", $sub->name, "\n";
   my $outfile=$sub->name . '.fas';
   $sub->printToFile($outfile);
   print "subsequence written to file: $outfile\n";
}

sub generateSub {
   chdir $workdir;
   my $guide=readGuide($guideFile);
   my @refseqFiles=("hcv1a_NC_004102.1.fas", "hcv1b_AJ238799.fas",
      "hcv2_AJ238799.fas", "hcv3_GU814263.fas", "hcv4_GU814265.fas",
      "hcv5_AF064490.fas", "hcv6_Y12083.fas");
   my %combinedbyns;
   #my @refseqFiles=glob("hcv*.fas");
   foreach my $f (@refseqFiles) {
      print "woking on $f ...\n";
      my $fasreader = new FastaReader($f);
      my $faseq = $fasreader->next;
      my $seqname=$faseq->name;
      print "guide for $seqname: ", join(' | ', @{$guide->{$seqname}}), "\n";
      push @{$combinedbyns{$subnames[0]}}, writeSub($faseq, $guide->{$seqname}->[0], $subnames[0]);
      push @{$combinedbyns{$subnames[1]}}, writeSub($faseq, $guide->{$seqname}->[1], $subnames[1]);
      push @{$combinedbyns{$subnames[2]}}, writeSub($faseq, $guide->{$seqname}->[2], $subnames[2]);
      push @{$combinedbyns{$subnames[3]}}, writeSub($faseq, combineRange($guide->{$seqname}->[1], $guide->{$seqname}->[2]),  $subnames[3]);
   }
}

sub combineRange {
   my $r1 = shift;
   my $r2 = shift; # number1-number2 in string format
   my @arr1 = split /-/, $r1;
   my @arr2 = split /-/, $r2;
   return $arr1[0] . '-' . $arr2[1];
}




__END__

=head1 DESCRIPTION

This is a quick and dirty method to generate the ference sequence for the DDL project.
May need more polishing.
