### perl ###

# to generate a table of request to be loaded into the database
# This program has to be run from the directory contain the mira runs

use FastaReader;
use Fastaseq;
#use Config::Simple;

my $rundir = '/d1/work/assembly';
my $projname; 
my $miraResultFile = "bestResult.tab";
my $resultDir = "results";

$outfile='miraresult.tab';

if (@ARGV < 1) { usage(); }
my $i = 0;
while ($ARGV[$i]) {
   if ($ARGV[$i] eq '-o') {
      $outfile = $ARGV[++$i];
   }
   else {
      $projname = $ARGV[$i];
   }
   ++$i;
}

# all work is done in the project subdirectory
chdir $projname;
my $miraResult = readResult($miraResultFile);
transformResult($miraResult);
chdir "-";
print "done\n";


#############################
sub transformResult {
   my $result = shift;
   my $seqreader = new FastaReader();
   open OU, ">$outfile" or die $!;
   my $h = $result->[0];
   splice @$h, 8, 2;
   push @$h, ("request", "sequence");
   print OU join("\t", @$h), "\n";
   for (my $i=1; $i < @$result; ++$i) {
      splice @{$result->[$i]}, 8, 2;
      #$result->[$i][0] = removeSuffix($result->[$i][0]);
      #print join("\t", @{$result->[$i]}), "\n";
      my $fileName;
      if ($result->[$i][1] == 1) {
         $fileName = $result->[$i][0] . "_c*.fas";
      }
      else {
         $fileName = $result->[$i][0] . '_' . fraction2string($result->[$i][1]) . "_*.fas";
      }
      my @files = glob("$resultDir/$fileName");
      if (@files == 1) {
         # should have only one file!
         #print "we got it: $files[0]\n";
      }
      elsif (@files > 1) {
         print "$fileName found ", join('  ', @files), "\n";
         die "found more than one sequence consensus\n";
      }
      else {
         die "$fileName found no sequence consensus, you may need to run runmira -c\n";
      }
      $seqreader->useFile($files[0]);
      my $seq = $seqreader->next;
      push @{$result->[$i]}, ($projname, $seq->sequence);
      print OU join("\t", @{$result->[$i]}), "\n";
   }
   system("pwd");
   print "result written to $outfile\n";
}

sub usage {
   print "buildrequest -d 'mira automatic assembly' -r 'Kemin Zhou' bar1\n";
   exit 1;
}

sub readResult {
   my $file = shift;
   open IN, "<$file" or die $!;
   my $line = <IN>; # header
   chomp $line;
   my @header = split /\t/, $line;
   if (@header != 11) {
      print scalar(@header), " columns not 11 as expected!\n";
      die;
   }

   my @results;
   push @results, \@header;
   while (<IN>) {
      chomp;
      my @row = split /\t/, $_;
      push @results, \@row;
   }
   return \@results;
}

sub removeSuffix {
   my $str = shift;
   if ($str =~ /(.+)\.fastq/) {
      return $1;
   }
   return $str;
}

sub fraction2string {
   my $frac = shift;
   if ($frac == 1) {
      return 1;
   }
   return int(sprintf("%.4f", $frac)*10000);
}

