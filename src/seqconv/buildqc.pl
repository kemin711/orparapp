### perl ###

use FastaReader;
use Fastaseq;
#use Config::Simple;

my $rundir = '/d1/work/assembly';
my $projname; 
my $qcFile = "QC.tab";
my $resultDir = "results";

$outfile='miraqc.tab';

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
my $miraqc = readQC($qcFile);
transformQC($miraqc);
chdir "-";
print "done\n";


#############################
sub transformQC {
   my $result = shift;
   open OU, ">$outfile" or die $!;
   my $h = $result->[0];
   print OU join("\t", @$h), "\n";
   for (my $i=1; $i < @$result; ++$i) {
      my $fileName;
      if ($result->[$i][1] == 1) {
         $fileName = $result->[$i][0] . "_c*.aln";
      }
      else {
         $fileName = $result->[$i][0] . '_' . fraction2string($result->[$i][1]) . "_*.aln";
      }
      my @files = glob("$resultDir/$fileName");
      my $correctFile;
      if (@files == 1) {
         $correctFile = $files[0];
         # should have only one file!
         #print "we got it: $files[0]\n";
      }
      elsif (@files == 2) {
         # one direction is wrong, used the other direction
         if ($result->[$i][3] eq '-') { # direction of ref reverse complement
            ($correctFile) = grep { /_rc1\.aln$/ } @files;
         }
         elsif ($result->[$i][4] eq '-') { # rc of assembly consensus
            ($correctFile) = grep { /_rc2\.aln$/ } @files;
         }
         else {
            die "Found two align file one is bad and the other is good\n",
               join(' | ', @files), "\n";
         }
         print "correct file: $correctFile out of the two: ", join(' | ', @files), "\n";
      }
      else {
         print "$fileName\n";
         die "did not find sequence alnfile\n";
      }
      my $alntxt = readAlnText($correctFile);
      push @{$result->[$i]}, $alntxt;
      print OU join("\t", @{$result->[$i]}), "\n";
   }
   system("pwd");
   print "result written to $outfile\n";
}

sub readAlnText {
   my $file = shift;
   open IN, "<$file" or die $!;
   my $buffer = '';
   while (<IN>) {
      chomp;
      $buffer .= ($_ . '\n');
   }
   close IN;
   return $buffer;
}


sub usage {
   print "buildrequest -d 'mira automatic assembly' -r 'Kemin Zhou' bar1\n";
   exit 1;
}

=head2 readQC

 header:
   input fraction ends direction_ref direction_ass alnlen numgap_ref numgap_ass gaplen_ref gaplen_ass identity begin_ref end_ref begin_ass end_ass alntext); # add other columns

=cut
sub readQC {
   my $file = shift;
   print "reading original QC from $file ...\n";
   open IN, "<$file" or die $!;
   my $line = <IN>; # header
   # preserve first 2 as natural key to find entry in result table
   chomp $line;
   my @header = split /\t/, $line;
   if (@header != 9) {
      print scalar(@header), " columns not 9 as expected!\n";
      die;
   }

   my @results;
   splice @header, 2, 4; # input_coverage ...
   splice @header, 3, 1; # remove direction
   pop @header; # remove last column
   push @header, qw(ends direction_ref direction_ass alnlen numgap_ref numgap_ass gaplen_ref gaplen_ass identity begin_ref end_ref begin_ass end_ass alntext); # add other columns
   push @results, \@header;
   while (<IN>) {
      chomp;
      my @row = split /\t/, $_;
      splice @row, 2, 4;
      splice @row, 3, 1;
      expandLast(\@row);
      push @results, \@row;
   }
   return \@results;
}

sub expandLast {
   my $r = shift;
   my $str = pop (@$r); # last column


   my @tmp = split / /, $str;
   my $refassPair = shift @tmp; # remove seqlen column
   my @pair = split /\|/, $refassPair;
   if (@pair != 2) {
      print "expanding last column $str ...\n";
      die "expecting pair of sequence names, $refassPair\n";
   }
   my @row = ();
   my ($refdir, $assdir);
   $refdir='+';
   $assdir='+';
   if ($pair[0] =~ /_rc$/ || $pair[0] =~ /rc$/) {
      $refdir = '-';
   }
   if ($pair[1] =~ /rc$/) {
      $assdir='-';
   }
   push @row, $refdir, $assdir;
   # discard sequence length, we should have them
   shift @tmp;

   # extract alnlen numgap_ref numgap_ass gaplen_ref gaplen_ass identity
   if ($tmp[0] =~ /alnlength=(\d+)/) {
      push @row, $1;
   }
   else { 
      print "you have ", $tmp[0], "\n";
      die "expecting alnlength\n"; 
   }

   if ($tmp[1] =~ /numgap=(\d+),(\d+)/) {
      push @row, $1, $2;
   }
   else { die "expecting numgap\n"; }

   if ($tmp[2] =~ /gaplength=(\d+),(\d+)/) {
      push @row, $1, $2;
   }
   else { die "expecting gaplength\n"; }

   if ($tmp[3] =~ /identity=([\d.]+)/) {
      push @row, $1;
   }
   else { die "expecting identity\n"; }

   if ($tmp[4] =~ /range1=(\d+)-(\d+)/) {
      push @row, $1, $2;
   }
   else { die "expecting range1\n"; }

   if ($tmp[5] =~ /range2=(\d+)-(\d+)/) {
      push @row, $1, $2;
   }
   else { die "expecting range2\n"; }

   push @$r, @row;
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

