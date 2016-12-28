### perl ###

# to generate a table of request to be loaded into the database
# This program has to be run from the directory contain the mira runs

use FastaReader;
use Fastaseq;
use Config::Simple;

my $rundir = '/d1/work/assembly';
my ($projname, $description, $requester, $outfile);

$requester='Kemin Zhou';
$description='Mira automatic assembly';
$outfile='miraproj.tab';

if (@ARGV < 1) { usage(); }
my $i = 0;
while ($ARGV[$i]) {
   if ($ARGV[$i] eq '-d') {
      $description = $ARGV[++$i];
   }
   elsif ($ARGV[$i] eq '-r') {
      $requester = $ARGV[++$i];
   }
   else {
      $projname = $ARGV[$i];
   }
   ++$i;
}

my $refseqFile = "$projname/$projname" . 'known.fas';
my $seqreader = new FastaReader($refseqFile);
my $faseq = $seqreader->next;

my $cutcfg = new Config::Simple("$projname/linearize.cfg");
my $headlen = length($cutcfg->param('head'));
my $taillen = length($cutcfg->param('tail'));
my $sequence = $faseq->sequence;
my @header = qw(name description head_len tail_len requester refseq);

# write the file
# name description head_len tail_len requester refseq
my $needHeader = 0;
if (! -f $outfile) {
   $needHeader = 1;
}
open OU, ">>$outfile" or die $!;
if ($needHeader) {
   print OU join("\t", @header), "\n";
}

print OU "$projname\t$description\t$headlen\t$taillen\t$requester\t$sequence\n";
print "output written to $outfile\ndone\n";

sub usage {
   print "buildrequest -d 'mira automatic assembly' -r 'Kemin Zhou' bar1\n";
   exit 1;
}
