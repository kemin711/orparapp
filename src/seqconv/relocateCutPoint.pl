### perl ###
use Fastaseq;
use FastaReader;

my ($infile, $cutpos);
my $i = 0;
while ($ARGV[$i]) {
   if ($ARGV[$i] eq "-i") {
      $infile = $ARGV[++$i];
   }
   elsif ($ARGV[$i] eq "-c") {
      $cutpos = $ARGV[++$i]; 
   }
   else {
      $infile = $ARGV[$i];
   }
   ++$i;
}


my $reader = new FastaReader($infile);
my $faseq = $reader->next;
my $newseq = $faseq->relocateCutPoint($cutpos);
my $outfile = $newseq->name . ".fas";
$newseq->printToFile($outfile);
print "new sequence written to $outfile\n";

__END__

=head1 NAME

relocateCutPoint - moving the cut point of circular DNA to a new position

=head1 SYNOPSIS

 relocateCutPoint -c 10 pasmid.fas

=head1 DESCRIPTION

This program relocated the cut point to a different position.

=head2 Options

 -i input plasmid fasta formated file
 -c new cut point, 1-based index in the old coordiante system
    The new version will start at new_cut_point

 The input file can also be given as the only argument.

=head1 AUTHOR

Kemin Zhou
