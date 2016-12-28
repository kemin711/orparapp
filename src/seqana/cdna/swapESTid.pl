#!/usr/bin/perl -w

use Bioseq;
use SeqIterator;

# to use EST id instead of the NCBI gi number 
# because some of the programs can guess the direction
# of the EST by the .x1 .y1 tag pairs
# >7413258 832002F09.y1 C. reinhardtii CC-125 nutrient replete, -S, -Fe,
# becomes
# >832002F09.y1

my $reader=SeqIterator->new($ARGV[0]) or die $!;
while (my $seq=$reader->next) {
    my $id=$reader->getId;
    my $newid=$id;
    my $title=$reader->title;
    if (defined $title && $title =~ /^([0-9A-Z]+\.[a-z]\d) (.+)/) {
        #print STDERR "got title $title\n";
        $newid=$1;
        $title= "gi_$id " . $2;
    }
    print ">$newid";
    if (defined $title) { print " $title"; }
    print "\n";
    printSeq($seq);
}
