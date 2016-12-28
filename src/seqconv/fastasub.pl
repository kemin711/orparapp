### perl ###
use Fastaseq;
use FastaReader;

my ($infile, $begin, $end);
$begin = 0;
$end = -1;
my $i = 0;
while ($ARGV[$i]) {
   if ($ARGV[$i] eq "-i") {
      $infile = $ARGV[++$i];
   }
   elsif ($ARGV[$i] eq "-b") {
      $begin = $ARGV[++$i]; 
   }
   elsif ($ARGV[$i] eq "-e") {
      $end = $ARGV[++$i]; 
   }
   elsif ($ARGV[$i] eq "--help") {
      usage();
   }
   else {
      $infile = $ARGV[$i];
   }
   ++$i;
}

my $reader = new FastaReader($infile);
my $faseq = $reader->next;
my $newseq = $faseq->subseq($begin, $end);
my $outfile = $newseq->name . ".fas";
$newseq->printToFile($outfile);
print "new sequence written to $outfile\n";

sub usage {
   print STDERR "fastasub -b 10 -e 99 inputfile.fas";
   exit 1
}


__END__

=head1 NAME

fastasub - take a subsequence from an fasta sequence file

=head1 SYNOPSIS

 fastasub -b 10 -e 99 inputfile.fas


=head1 SEE ALSO

seqed in the ../seqana directory was an earlier program.
The seqed program is old and needs to be polished


