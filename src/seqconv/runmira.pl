### perl ###
my $HOSTNAME;
BEGIN {
   $HOSTNAME=$ENV{HOSTNAME};
   print "hostname: $HOSTNAME\n";
   if ($HOSTNAME eq 'foo') {
      $ENV{ORACLE_HOME} = "/u01/home/oracle/product/11.2.0/dbhome_1";
      my @extraPath = ("/apps/sys/perl-5.18.1/lib/site_perl/5.18.1", 
            "/apps/sys/perl-5.18.1/lib/site_perl/5.18.1/x86_64-linux/auto/DBD",
            "/apps/sys/perl-5.18.1/lib/site_perl/5.18.1/x86_64-linux",
            "/home/zhouke/foo/perlmod",
            "/home/zhouke/foo/perlmod/lib/site_perl/5.18.1");
      foreach my $p (@extraPath) {
         if (! grep { $_ eq $p } @INC) {
            print "adding $p to perl5lib\n";
            push @INC, $p;
         }
         else {
            print "$p already in perl5lib\n";
         }
      }
      $ENV{LD_LIBRARY_PATH}="/apps/sys/gcc-4.9.0/lib64:/apps/sys/lib:/apps/sys/perl-5.18.1/lib/site_perl/5.18.1/x86_64-linux/auto/DBD/Oracle:/apps/sys/perl-5.18.1/lib/site_perl/5.18.1/x86_64-linux/auto/DBI";
   }
   print "Perl5 Environment: ", join(" | ", @INC), "\n";
}

use strict;
use Log::Log4perl;
use FastaReader;
use Fastaseq;
use MiraResult;
use Config::Simple;
use File::Basename;
use CutDetector;
use Pod::Usage;
use DBI;

#my $miraBinary = "/usr/local/bioinfo/bin/mira"; # this is for europa
#my $binaryDir = "/usr/local/bioinfo/bin";
my $tmpdir = '/d1/tmp';

my $miraBinary;
my $binaryDir; #= "/usr/local/bioinfo/bin";
if ($HOSTNAME eq 'foo') {
   $miraBinary = "/home/zhouke/foo/download/mira_4.0rc3_linux-gnu_x86_64_static/bin/mira";
   $tmpdir = '/scratch/mirarun';
   $binaryDir="/home/zhouke/foo/bin";
   if (index($ENV{PATH}, $binaryDir) == -1) {
      $ENV{PATH} = $ENV{PATH} . ":$binaryDir";
   }
}
elsif ($HOSTNAME eq 'bar') {
   $miraBinary = "/usr/local/bioinfo/mira_4.0rc3_linux-gnu_x86_64_static/bin/mira";
   $binaryDir = "/usr/local/orpara/bin";
   #$ENV{ORACLE_HOME} = "/data/kzoracle/product/12.1.0/dbhome_1";
   #$ENV{LD_LIBRARY_PATH}="/data/kzoracle/product/12.1.0/dbhome_1/lib:/usr/local/lib:/usr/local/lib64:/usr/lib64:/usr/local/lib64/perl5/auto/DBD/Oracle";
   #$ENV{PERL5LIB}="/usr/local/lib64/perl5/auto/DBD:/usr/lib/perl5/site_perl:/home/zhouke/perlmod:.";
}
print "mira binary: $miraBinary\n";

# for kraken we will use a different location.
# right now this pipeline is based on the 4.0.1rc3
# rc5 and stable release 4.0.1 have quite different behaviors
# about 1-2 weeks testing and code change is needed to use the
# latest version of mira

#my $STANDARD_START_BASE=4000000; # fixed amount
#my $STANDARD_PLASMID_LENGTH=9000; # the fixed plasmid length for START_BASE
my $HIGH_COVERAGE = 1700;
my $LOW_COVERAGE = 100;
my $ASCENT=3;
# this is the difference between the neighboring runs
# used by findHighestFraction. It should be ralative to
# the range, better 1/10 of the range rather than a
# fixed number.
my $STOP_CONDITION=0.01;
my $STEP = 0.1;
# for QC, when display the alignment, flip the consensus
# if the consensus is the RC of the reference.
# If set to false, it will flip the reference.
# Flippping the consensus makes comparing the difference
# easy, but makes looking at the assembly through GAP5 
# a tiny bit harder (You have to do length - coordinates)
my $FLIP_CONSENSUS = 0;

my $plasmidLength=9000;
#my $startWith = 4000000; # start with this if plasmi is 9kb
# this should be a fixed amount, only length and average length
# should vary

# the lower starting point, guess, that gives a good assembly
# but not too much. Starting point for the lower end search.
# Too high may be beyong the limit of Mira can handle, too low
# will also not be enough to generate a continus contig.
#
#my $inputseq = "bar2.fastq";
my $action=0;
my $i=0;
my $inputseq;

readConfig();
while ($ARGV[$i]) {
   #if ($ARGV[$i] eq '-s') { $STANDARD_START_BASE = $ARGV[++$i]; }
   if ($ARGV[$i] eq '-l') { $plasmidLength = $ARGV[++$i]; }
   elsif ($ARGV[$i] eq '-H') { $HIGH_COVERAGE = $ARGV[++$i]; }
   elsif ($ARGV[$i] eq '-L') { $LOW_COVERAGE = $ARGV[++$i]; }
   elsif ($ARGV[$i] eq '--flip-consensus') { $FLIP_CONSENSUS = 1; }
   elsif ($ARGV[$i] eq '--flip-reference') { $FLIP_CONSENSUS = 0; }
   elsif ($ARGV[$i] eq '-A') { $ASCENT = $ARGV[++$i]; }
   elsif ($ARGV[$i] eq '--clean') { clean(); }
   elsif ($ARGV[$i] eq '-c') { 
      $action = 1;
   }
   elsif ($ARGV[$i] eq '--help') { pod2usage(1); }
   elsif ($ARGV[$i] =~ /^-/) {
      die "option $ARGV[$i] not able to process\n";
   }
   else {
      $inputseq = $ARGV[$i];
   }
   ++$i;
}

Log::Log4perl->init("/home/zhouke/src/proj/seqconv/log.conf");
my $log = Log::Log4perl->get_logger("runmira");
$log->info("===========================");
$log->info("| Starting a new mira run |");
$log->info("===========================");
$log->info("plasmid length: $plasmidLength nt, ASCENT: $ASCENT");

my $cutcfg = new Config::Simple("linearize.cfg");
if ($action == 1) {
   checkResult();
   exit 0;
}

if (!$inputseq) {
   die "You need an input sequence\n";
}

my @trimParameters=( [4,20], [4,19], [5,20], [5,19], [6, 20], [6,19] );
# cut site should be stored in linearize.cfg file in each directory
#my $cut_site="gacggatcgggagatctcccgatcccctatggtgcactctcagtacaatctgct";

my $uncutResults;
open OU1, ">bestResult.txt" or die $!;
open OU2, ">bestResult.tab" or die $!;
eval {
   my ($numGoodContig, $numCut, $cutResults) = assembleCut();
   if (@$cutResults > 0) {
      showResult($cutResults, \*OU1);
      writeResult($cutResults, \*OU2, 1, 1);
   }
   # need to do uncut assembly
   if ($numCut > 1 || $numCut > $numGoodContig) {
      $uncutResults = assembleUncut();
      showResult($uncutResults, \*OU1);
      writeResult($uncutResults, \*OU2, 0, 0);
   }
};
if ($@) {
   $log->fatal("the mira pipeline failed with error: $@");
   stampFail();
}
else {
   $log->info("best rbestResult.txtesult written to bestResult.txt");
   if (updateStatusDone()) {
      stampSuccess();
   }
   else {
      stampFail();
   }
}
close OU1;
close OU2;

################### sub routines #################################

sub stampSuccess {
   if (-e "runmira.fail") {
      unlink "runmira.fail";
   }
   system("touch runmira.done");
   $log->info("runmira done");
}

sub stampFail {
   if (-e "runmira.done") {
      unlink "runmira.done";
   }
   system("touch runmira.fail");
   $log->fatal("runmira fail");
}

sub clean {
   #system("rm bestResult.*"); save the result file, very samll
   #and useful
   system("rm -rf *cut* *trim*");
   system("rm *.done");
   system("rm *.log");
   exit(0);
}

sub readConfig {
   if (-f 'runmira.cfg') {
      my $config = new Config::Simple("runmira.cfg");
      if ($config->param('input')) {
         $inputseq = $config->param('input');
      }
      if ($config->param('length')) {
         $plasmidLength = $config->param('length');
      }
      if ($config->param('high')) {
         $HIGH_COVERAGE = $config->param('high');
         $LOW_COVERAGE = $config->param('low');
      }
      if ($config->param('ascent')) {
         $ASCENT=$config->param('ascent');
      }
      return 1;
   }
   else {
      warn "There is no config file runmira.cfg\n";
      return 0;
   }
}

=head2 checkResult

Display some basic info.
Then align result to reference.
Will also fetch raw sequence to the result directory.

If run successful, it should load result to database.
At this point, I am not sure it should be done in
this program or miraload.

=cut
sub checkResult {
   my $head = uc($cutcfg->param('head'));
   my $tail = uc($cutcfg->param('tail'));
   open CK, ">QC.tab" or die $!;
   print CK "input\tfraction\tinput_coverage\tlength\tquality\tcoverage\tends\tdirection\tcomparison_to_reference\n";

   print $log->debug("tail-head $tail | $head");
   if (-d "results") {
      system("rm -f results/*");
   }
   else {
      mkdir "results";
   }
   my $references = prepareReference($inputseq, "results");

   # HEADER for bestResult.tab file
   # 0     1         2               3          4         5    6    7    8  
   #input frac input_coverage length_largest quality coverage head tail head_rc  
   #   9       10
   # tail_rc  rawcut
   my $results = readResult(); # reading from bestResult.tab
   $log->info(scalar(@$results) . " results to check");
   shift @$results; # get rid of header
   foreach my $r (@$results) {
      $r->[1] =~ s/0+$//;
      $log->debug("working on [" . $r->[0] . " " . $r->[1] . "]");
      my $good = 1;
      for (my $j = 6; $j <= 9; ++$j) { # upper case sequence
         $r->[$j] = uc($r->[$j])
      }
      $r->[2]=sprintf("%i", $r->[2]); # format input coverage

      my @tmpArr = @{$r}[0..5];
      #$tmpArr[0] = removeSuffix($tmpArr[0]); # already removed, no need
      print CK join("\t", @tmpArr), "\t";
      # check the ends with cut point 
      # This is only meaningful for cut assembly
      my $rc = 0;
      my $endMatch = 0;
      if ($r->[1] == 0 || !$r->[6] || !$r->[7]) {
         $log->debug("Assembly bad no QC info");
         $good = 0;
      }
      if ($r->[10] && $good) { # raw sequence cut
         if ($head eq $r->[6] && $tail eq $r->[7]) {
            $endMatch = 2;
         }
         elsif ($tail eq $r->[8] && $head eq $r->[9]) {
            $endMatch = 2;
            $rc = 1;
         }
         elsif ($head eq $r->[6] || $tail eq $r->[7]) {
            $endMatch = 1;
         }
         elsif ($tail eq $r->[8] || $head eq $r->[9]) {
            $endMatch = 1;
            $rc = 1;
         }
         else {
            $rc = -1; # unknown
         }
      }

      print CK "$endMatch\t";
      if ($good) {
         my ($refseqFileName, $extract_rc);
         if ($FLIP_CONSENSUS) {
            $refseqFileName = $references->[0];
            $extract_rc = $rc;
         }
         else {
            if ($rc != -1) {
               $refseqFileName = $references->[$rc];
            }
            else {
               $refseqFileName = $references->[0];
            }
            $extract_rc = 0;
         }
         my $sfn = extractSequence($r->[0], $r->[1], "results", $extract_rc);

         my ($alnrc, $alninfo);
         if ($r->[10]) {
            ($alnrc, $alninfo) = alignToRefCut($refseqFileName, $sfn);
            if ($alnrc == 1) {
               $rc = !$rc;
            }
         }
         else {
            # the direction of the uncut assembly is not known
            ($rc, $alninfo) = alignToRefUncut($refseqFileName, $sfn);
         }
         print CK directionToString($rc), "\t$alninfo";

      }
      print CK "\n";
   }
   close CK;
   print "result written to QC.tab\n";
   exit 0;
}

sub directionToString {
   my $val = shift;
   if ($val == 0) { return "ref"; }
   elsif ($val == 1) { return "rc_ref"; }
   else { return "unknown"; }
}

# the logic here is not very well thought out
# May need factoring in the future.
sub alignToRefUncut {
   my $refname = shift;
   my $fasfile = shift;

   chdir "results";
   my ($rc, $info) = alignTwoSeq($refname, $fasfile);
   # find new cut point for reference sequence
   my ($rb1, $re1) = extractRange($info, 1);
   # read reference sequence into object

   my ($refCutFile, $dummy);
   if ($FLIP_CONSENSUS || !$rc) {
      $refCutFile = recutFile($refname, $rb1, $re1);
      ($dummy, $info) = alignTwoSeq($refCutFile, $fasfile, $rc);
   }
   else { # flip Reference and RC == True
      my $rcRefname = getRCReferenceFileName($refname);
      $refCutFile = recutFile($rcRefname, $rb1, $re1);
      ($dummy, $info) = alignTwoSeq($refCutFile, $fasfile);
   }

   chdir "..";
   return ($rc, $info);
}

sub recutFile {
   my $file = shift;
   my $b = shift;
   my $e = shift;

   my $newfas;
   my $reader = new FastaReader($file);
   my $fas = $reader->next;
   if ($b == 1) {
      # range1=1-5198 range2=3943-9914
      #     1 
      #     =========> Ref
      # =========>     Ass
      $newfas = $fas->relocateCutPoint($e+1);
   }
   else {
      # ==============>
      #     ============>
      $newfas = $fas->relocateCutPoint($b);
   }
   my $newfile = $newfas->name . ".fas";
   if (!-f $newfile) {
      $newfas->printToFile($newfile);
   }
   return $newfile;
}

sub getRCReferenceFileName {
   my $rf = shift;
   if ($rf =~ /known\.fas/) {
      $rf =~ s/known\.fas/knownrc\.fas/;
   }
   elsif ($rf =~ /knowrc\.fas/) {
      $rf =~ s/knownrc\.fas/known\.fas/;
   }
   return $rf;
}

sub extractRange {
   my $info = shift;
   my $number = shift;

   my ($rb, $re);
   if ($info =~ /range$number=(\d+)-(\d+)/) {
      $rb = $1;
      $re = $2;
   }
   else { die "alninfo in wrong format: $info\n"; }
   return ($rb, $re);
}

sub extractRangeBoth {
   my $info = shift;

   my ($rb1, $re1, $rb2, $re2);
   if ($info =~ /range1=(\d+)-(\d+)/) {
      $rb1 = $1;
      $re1 = $2;
   }
   else { die "alninfo in wrong format: $info\n"; }
   if ($info =~ /range2=(\d+)-(\d+)/) {
      $rb2 = $1;
      $re2 = $2;
   }
   else { die "alninfo in wrong format: $info\n"; }
   return ($rb1, $re1, $rb2, $re2);
}

=head2 alignToRef

Align one assembly to reference sequence.
Match REF to Assembly

=cut
sub alignToRefCut {
   my $refname = shift;
   my $fasfile = shift;

   chdir "results";
   my ($rc, $alninfo) = alignTwoSeq($refname, $fasfile);
   chdir "..";
   return ($rc, $alninfo);
}

sub alignTwoSeq {
   my $fileRef = shift;
   my $fileAss = shift;
   my $rc = shift;

   my $fasname = basename($fileAss);
   my $stem = getFileStem($fasname);
   my $outfile = $stem . ".aln";
   my $alninfo;

   if (!defined $rc) { $rc = 0; }
   if ($rc) {
      if ($FLIP_CONSENSUS) {
         $outfile = $stem . "_rc2.aln";
         $alninfo = `alnlocal -o $outfile -r 2 $fileRef $fasname | head -1`;
      }
      else {
         $outfile = $stem . "_rc1.aln";
         $alninfo = `alnlocal -o $outfile -r 1 $fileRef $fasname | head -1`;
      }
   }
   else {
      $alninfo = `alnlocal -o $outfile $fileRef $fasname | head -1`;
      if ($alninfo =~ /identity=([0-9.]+)/) {
         my $identity = $1;
         if ($identity < 0.82) {
            $log->debug("bad alignment: $alninfo, try the other strand");
            if ($FLIP_CONSENSUS) {
               $outfile = $stem . "_rc2.aln";
               $alninfo = `alnlocal -o $outfile -r 2 $fileRef $fasname | head -1`;
            }
            else {
               $outfile = $stem . "_rc1.aln";
               $alninfo = `alnlocal -o $outfile -r 1 $fileRef $fasname | head -1`;
            }
            $rc = 1;
         }
      }
      else {
         die "something is wrong with the alnlocal program stdout output: $alninfo\n";
      }
   }
   chomp $alninfo;
   #if ($alninfo =~ s/^.+?\s//) { # shorten remove the first few 
   #}
   #else {
   #   die "alnlocal infor line wrong $alninfo\n";
   #}
   return ($rc, $alninfo);
}

sub copyReferenceRC {
   my $refpath = shift;
   my $dir = shift;
   my $outname = basename($refpath);
   $outname =~ s/\.fas/rc.fas/;
   if (-f $outname) {
      return $outname;
   }

   # this should be used in the working directory
   my $reader = new FastaReader($refpath);
   my $reffas = $reader->next;
   $reffas->revcomp;
   $reffas->printToFile("$dir/$outname");
   print "rc reference file $outname written to directory $dir\n";
   return $outname;
}

sub prepareReference {
   my $inputFastqFile = shift;
   my $dir = shift;

   my $refs = buildReferenceName($inputFastqFile);
   copyReferencesTo($refs, $dir);
   return $refs;
}

sub buildReferenceName {
   my $inputFastq = shift;
   my $ref = getFileStem($inputFastq) . "known.fas";
   my $rcRef = getFileStem($inputFastq) . "knownrc.fas";
   return [$ref, $rcRef];
}

sub copyReferencesTo {
   my $references = shift;
   my $dir = shift;

   if (-e $references->[0]) {
      system("cp", $references->[0], $dir);
      system("chmod +w ", $dir . "/" . $references->[0]);
   }
   else {
      $log->fatal("unix cp command source file " . $references->[0] . " missing\n");
      die "Cannot find source file for unix cp command ", $references->[0], "\n";
   }
   my $reader = new FastaReader($references->[0]);
   my $reffas = $reader->next;
   $reffas->revcomp;
   my $rcFile = $dir . "/" . $references->[1];
   $reffas->printToFile($rcFile);
   print "rc reference file ", $references->[1], " written to directory $dir\n";
}

=head extractSequence 

 parameter($fileName, $fraction, $directory, $rc)
  $fileName is the name of the input fastq file for doing the assembly
  $fraction is the fraction used for the assembly.
  $directory where to put the sequence file into.
  $rc whether to reverse-complement the sequence or not.

=cut
sub extractSequence {
   my $file = shift; # file name, no longer has the .fastq suffix
   my $frac = shift; # best fraction
   my $dir = shift;
   my $rc = shift; # 1 for revcomp, 0 for forward, -1 for unknown

   $log->debug("extracting consesus sequence from file $file at fraction $frac into directory $dir");

   my $partialFileName;
   if ($frac == 1) {
      $log->debug("100% not partial file");
      $partialFileName = $file;
   }
   else {
      $log->debug("$frac partial file");
      $partialFileName = makePartialFileName($file, $frac);
   }

   my $stem = getFileStem($partialFileName);
   my $miraResult = MiraResult->new($stem);
   $miraResult->evaluateResult;
   my $seq = $miraResult->getLargestContigSequence;
   if ($rc == 1) {
      $seq->revcomp;
   }
   my $seqFileName = $dir . "/" . $seq->name . ".fas";
   $seq->printToFile($seqFileName, 80);
   return $seqFileName;
}

=head2 removeSuffix

Only removed the fastq suffix

=cut
sub removeSuffix {
   my $in = shift;
   if ($in =~ /(.+)\.fastq/) {
      return $1;
   }
   else { return $in; }
}
=head assembleCut

   return ($foundOneContig, $naturalCut, \@cutResult);

   cutResult: [$cutfile, [ @goodResults ] ]

=cut
sub assembleCut {
   my @cutResult = ();
   my $foundOneContig = 0;
   my $naturalCut = 0;
   my $cutfile = linearizeInput($inputseq);
   my ($turnFrac, $mres, $CPF) = getTurningFraction($cutfile);
   my ($allResults);
   if ($turnFrac == 0) {
      ++$naturalCut;
      $log->info("Could not find a fraction to assembly cut sequence into one contig");
   }
   elsif ($turnFrac == -1) {
      return (0, 1, \@cutResult);
   }
   else {
      $log->info("Turning fraction from cut input: $turnFrac, with CPF $CPF");
      $foundOneContig = 1;
      my @goodResults;
      #                    frac      result      coverage
      #if ($mres) { # has good turning point
      #   push @goodResults, [$turnFrac, $mres, $turnFrac*$CPF]; # top one
      #}
      $allResults = collectAllGoodMiraResults($cutfile, $CPF);
      push @goodResults, @$allResults;
      push @cutResult, [ $cutfile, [ @goodResults ] ];
   }

   # try trim and cut
   my $trimcutFiles = trimAndLinearize($inputseq);
   $log->info("trying trim-and-cut: " . join(", ", @$trimcutFiles) . " ...");
   foreach my $tf (@$trimcutFiles) {
      ($turnFrac, $mres, $CPF) = getTurningFraction($tf);
      if ($turnFrac > 0) {
         ++$foundOneContig;
         my @tmpGood=();
         #if ($mres) {
         #   push @tmpGood, [$turnFrac, $mres, $turnFrac*$CPF];
         #}
         $allResults = collectAllGoodMiraResults($tf, $CPF);
         push @tmpGood, @$allResults;
         push @cutResult, [ $tf, [ @tmpGood ] ];
      }
      elsif ($turnFrac == 0) {
         $log->warn("did not find turning fraction for $tf");
         ++$naturalCut;
      }
      else {
         $log->warn("did not find turning fraction for $tf");
      }
   }
   $log->info("$foundOneContig cut assemblies produced one large contig, $naturalCut got natural cuts");
   return ($foundOneContig, $naturalCut, \@cutResult);
}

=head2 assemblUncut

return true if need to do cut assembly.

=cut
sub assembleUncut {
   my @uncutResult = ();
   # 1. untrimmed
   my ($turnFrac, $mres, $CPF) = getTurningFraction($inputseq);
   my @good;
   if ($turnFrac > 0) {
      $log->info("Best fraction from raw uncut input: " . sprintf("%.4f", $turnFrac));
      #push @good, [$turnFrac, $mres, $turnFrac*$CPF];
      my $lowerResults = collectAllGoodMiraResults($inputseq, $CPF);
      push @good, @$lowerResults;
      push @uncutResult, [ $inputseq, [ @good ] ];
   }

   my $trimmed = trimsequences($inputseq);
   $log->debug("trimmed file: " . join(", ", @$trimmed) . "  need to be processed");
   foreach my $tf (@$trimmed) {
      my ($turnFrac, $mres, $CPF) = getTurningFraction($tf);
      if ($turnFrac > 0) {
         @good = ();
         #push @good, [$turnFrac, $mres, $turnFrac*$CPF];
         my $lowerResults = collectAllGoodMiraResults($tf, $CPF);
         push @good, @$lowerResults;
         push @uncutResult, [ $tf, [ @good ] ];
      }
   }
   return \@uncutResult;
}


=head2 collectAllLower

collect all the lower mira results

CPF = $HCoverage/$HFraction. A unit value to 
simplify calculation.

I am replacing it with collectAllGoodMiraResults.

=cut
sub collectAllLowerMiraResult {
   my $infile = shift;
   my $fraction = shift; # top one 
   my $CPF = shift;

   my @LLL = ();
   if ($fraction <= 0) {
      $log->debug("Nothing to do for $infile fraction at $fraction");
      return \@LLL;
   }

   #my $millage = int(10000*$fraction);
   # basis point 1/10000
   my $bips = int(10000*$fraction);
   my ($stem, $path, $suffix) = fileparse($infile, ".fastq");
   my $pattern = $stem . "_[1-9]*.fastq"; # for shell
   my @lowerFiles = `ls -r $pattern`;
   if (!scalar(@lowerFiles)) {
      $log->warn("There is no lower fractions run");
      return \@LLL;
   }
   chomp @lowerFiles;
   $pattern = $stem . "_([0-9]{1,5}).fastq"; # for perl
   foreach my $f (@lowerFiles) {
      if ($f =~ /$pattern/) {
         my $value = $1;
         if ($value < $bips) {
            $stem = getFileStem($f);
            my $mr = new MiraResult($stem);
            $mr->evaluateResult;
            if ($mr->getNumberOfGoodContigs == 1) {
               my $frac = $value/10000;
               #push @$good, [$frac, $mr, $frac*$CPF];
               push @LLL, [$frac, $mr, $frac*$CPF];
            }
            else {
               $log->warn("$f has no single good contig");
            }
         }
      }
      else {
         die "unexpected file pattern $f\n";
      }
   }
   @LLL = sort { $b->[0] <=> $a->[0] } @LLL;
   return \@LLL;
}

=head2 collectAllGoodMiraResults

 result array: fraction, MIRA_Result, fraction*CPF = coverage

=cut
sub collectAllGoodMiraResults {
   my $infile = shift;
   my $CPF = shift;

   my @LLL = ();

   my ($stem, $path, $suffix) = fileparse($infile, ".fastq");
   my $pattern = $stem . "_[1-9]*_assembly"; # for shell
   my @assemblies = glob($pattern);
   my $allUsed = $stem . "_assembly";
   if (-d $allUsed) {
      push @assemblies, $allUsed;
   }
   if (!scalar(@assemblies)) {
      $log->warn("failed to find any fastq file with pattern $pattern");
      return \@LLL;
   }
   $log->info(scalar(@assemblies) . " assemblies for stem $stem: " . join('  ', @assemblies));
   $pattern = $stem . "_([0-9]{1,5})_assembly"; # for perl
   foreach my $f (@assemblies) {
      if ($f =~ /$pattern/ || $f eq $allUsed) {
         my $value = $1; # the fraction in the file name
         my $projname = $stem . "_$value";
         if ($f eq $allUsed) {
            $projname = $stem;
            $value = 10000;
         }
         $log->debug("looking for mira result for mira project: $projname");
         eval {
            my $mr = new MiraResult($projname);
            $mr->evaluateResult;
            if ($mr->getNumberOfGoodContigs == 1) {
               my $frac = $value/10000;
               $log->debug("computed fraction: $frac");
               push @LLL, [$frac, $mr, $frac*$CPF];
            }
            else {
               $log->info("$f has no single good contig, ignored.");
            }
         };
         if ($@) {
            $log->error($@ . "mira run for $projname is ignored because either not done or failed somehow.");
         }
      }
      else {
         die "unexpected file pattern $f against $pattern\n";
      }
   }
   @LLL = sort { $b->[0] <=> $a->[0] } @LLL;
   $log->info(scalar(@LLL) . " good mira runs");
   return \@LLL;
}

sub showResult {
   my $bestf = shift;
   my $fh = shift;
   if (!$fh) { $fh = \*STDOUT; }
   
   foreach my $e (@$bestf) {
      my $good = $e->[1]; # fraction, MiraResult, coverage
      print $fh $e->[0], "\n";
      foreach my $g (@$good) {
         print $fh "    fraction: ", sprintf("%.5f", $g->[0]), " input coverage: ",
               sprintf("%i", $g->[2]), "\n    ";
         if ($g->[1]) { # has MiraResult
            $g->[1]->showBasicInformation($fh);
            if ($g->[0] != 0) {
               $g->[1]->showBreakPoint($fh);
            }
            else {
               print $fh "   poor assemble results not showing break point\n";
            }
         }
      }
      print $fh "\n";
   }
}

=head2 writeResult

Tabular format for the result summary.
each combination of parameters

 parameters(goodResult, FH_output, raw_cut_or_not, withHeaderOrNot)

=cut
sub writeResult {
   my $bestf = shift;
   my $fh = shift;
   my $rawcut= shift;
   my $withHeader = shift;

   if (!$fh) { $fh = \*STDOUT; }
   if (!defined $rawcut) { $rawcut = 0; }

   my $cutSiteLength = length($cutcfg->param('head'));

   $log->debug("writing best result from $bestf to tabular format");
  
   if ($withHeader) {
      # header, assembly changed to input 
      # Now writing the file extension name for input
      # the head/tail is simply the end cutSiteLength, if there is one or a few
      # bases extension, then it will be quite different from the original cut
      # site head tail
      print $fh "input\tfrac\tinput_coverage\tlength_largest\tquality\tcoverage\thead\ttail\thead_rc\ttail_rc\trawcut\n";
   }

   foreach my $e (@$bestf) { # one row of result [cutfile, [ @goodresults ]]
      foreach my $g (@{$e->[1]}) { #good results
         # $g structure: [fraction, MIRA_Result, fraction*CPF]
         # the number of digits has to be 2 longer than the file making
         # otherwise, there might be a problem with certain numbers.
         # The behavior of sprintf("%.4f", 0.20734517) is 0.2073
         # but sprintf("%.4f", 0.20735) is 0.2074.
         # So when you store sprintf("%.5f", 0.20734517), then read and 
         # recast to sprintf("%.4f", ) the value differe by the last digit!
         $log->debug("writing good result ", $e->[0] . " frac: " . $g->[0]);
         print $fh removeSuffix($e->[0]), "\t", sprintf("%.7f", $g->[0]), "\t", 
            sprintf("%.2f", $g->[2]), "\t";
         if ($g->[1]) { # has good assembly result
            print $fh $g->[1]->getLongestContig, "\t", $g->[1]->getQuality, "\t", 
                     $g->[1]->getCoverage, "\t"; 
            if ($g->[0] != 0) {
               print $fh join("\t", @{$g->[1]->getEndSequences($cutSiteLength)});
            }
            else {
               print $fh "\t\t\t\t";
            }
            print $fh "\t$rawcut";
         }
         else {
            $log->debug($e->[0] . " " . $g->[0], " has no good mira result");
            print $fh "\t\t\t\t\t\t\t\t$rawcut";
         }
         print $fh "\n";
      }
   }
}

=head2 readResult

Read the result from bestResult.tab file as a table
array of array data structure.

=cut
sub readResult {
   my $bestf = shift;
   #$log->info("Reading best results from bestResult.tab ...");

   my @best = ();
   open BE, "<bestResult.tab" or die "there is no bestResult.txt file\n";
   my $count = 0;
   my $line = <BE>;
   while ($line) {
      ++$count;
      #print $line;
      chomp $line;
      my @row = split /\t/, $line;
      push @best, [@row];
      $line = <BE>;
   }
   $log->info("read $count assemblies");
   return \@best;
}


=head2 getTurningFraction

 parameter($inputFile)

 Find the highest fraction of the input that still give one 
 large contig. This may not be the fraction given the best
 alignment.

 @return (fractionTurningPoint, MiraResult, CPF)
   The turningFraction will be 1 if not sufficient input.

 CPF = HCov/HFrac unit value
 coverage per fraction  CPF * fraction = Coverage

=cut
sub getTurningFraction {
   my $inputFile = shift; # fastq file containing all raw reads

   my $turnFrac = 0;
   my $result = undef;
   my ($mean, $count, $total) = evaluateFastq($inputFile);
   my ($topstatus, $lowerstatus, $topresult, $lowerresult);

   $log->info("$inputFile has $count fastq sequences, mean $mean, total $total bases");
   # HIGH_COVERAGE is just a hint, this program should discover the proper one
   my ($HFrac, $HCov) = coverageToFraction($HIGH_COVERAGE, $total);
   my $LCov = lengthAdjust($mean)*$LOW_COVERAGE;
   my $LFrac = $LCov*$plasmidLength/$total;
   if ($HCov < 350) {
      $log->warn("Insufficient coverage: " . sprintf("%i", $HCov) . " at $HFrac");
   }
   if ($LCov >= 0.9*$HCov) {
      $log->warn("Low coverage: " . $LCov . " too high. Cannot find a turning point");
      $LFrac = 0.8*$HFrac;
      $LCov = $LFrac * $total/$plasmidLength;
      $log->warn("Lowered lower coverage to " . $LCov . " Lower frac to $LFrac");
      #die "Lower coverage should not be higher than 0.9*HCov\n";
   }
   my $StartLFrac = $LFrac; # save away in case we need to change it later
   $log->info("Starting with fraction [" . sprintf("%.4f", $LFrac) . " - "
         . sprintf("%.4f", $HFrac) . "] and coverage: " . $LCov . " - " . $HCov);

   ($topstatus, $topresult) = runAndEvaluateMiraFraction($inputFile, $HFrac);
   my $delta = ($HFrac - $LFrac)/5;
   while ($topstatus == 1 && $HFrac < 1) { # adjust higher if not high enough
      $HFrac += $delta;
      if ($HFrac > 1) { $HFrac = 1; }
      $HCov = $HFrac*$total/$plasmidLength;
      $log->info("TOP COVERAGE: $HIGH_COVERAGE configured too low, increasing it to $HCov");
      ($topstatus, $topresult) = runAndEvaluateMiraFraction($inputFile, $HFrac);
   }
   if ($topstatus == 1) { # done, 
      if ($HFrac == 1) { # used all input, need to get a few more data points
         # and should try a few below to not miss the best assembly
         $log->info("Since reached good assembly with 100%, trying a few lower ones");
         my $lowres;
         for (my $f = 0.9; $f > 0.4; $f -= 0.1) {
            ($lowerstatus, $lowres) = runAndEvaluateMiraFraction($inputFile, $f);
            if ($lowerstatus != 1) {
               $log->debug("at $f fraction, no longer producing good assembly, stop now");
               last; # stop doing it if not producing good results
            }
         }
         $turnFrac = 1;
      }
      else {
         $log->fatal("Top coverage selected gave one good contig, consider increse your top coverage");
         die "Terminated due to TOP coverage too low\n";
      }
   }
   else { # start from the lowest fraction, 
      my $detector = new CutDetector;
      ($lowerstatus, $lowerresult) = runAndEvaluateMiraFraction($inputFile, $LFrac);
      $detector->absorb($lowerstatus);
      my $lastRun = 1;
      while ($LFrac < 0.9*$HFrac && $lowerstatus != 1) {
         $LFrac += $delta;
         $detector->absorb($lowerstatus);
         $log->debug("Increasing the lower fraction to $LFrac");
         if ($LFrac >= 0.99*$HFrac) {
            $lastRun = 0;
            last;
         }
         ($lowerstatus, $lowerresult) = runAndEvaluateMiraFraction($inputFile, $LFrac);
      }
      if ($lastRun) {
         $detector->absorb($lowerstatus);
      }
      $log->debug("detector state: " . $detector->state);
      if ($detector->hasCut) {
         $log->info("raw reads has natural cut");
         $turnFrac = 0;
      }
      elsif ($detector->hasNoCut) {
         $log->info("raw reads has no natural cut");
         $STOP_CONDITION = ($HFrac - $LFrac)/100;
         $STEP = ($HFrac - $LFrac)/10;
         $log->info("search turning frac in [$LFrac - $HFrac] ...\n");
         ($turnFrac, $result) = findHighestMiraFraction($inputFile, $LFrac, $HFrac, $lowerresult, $topresult);
      }
      elsif ($detector->hasError) {
         $log->info("detector in error state");
         die "state ", $detector->state, " is in error\n";
      }
      else {
         $log->info("don't know");
      }
      # special condition where we don't have enough input
      if ($detector->hasCut || $detector->hasError) {
         $log->warn("Could not find lower starting point from $StartLFrac to $HFrac, there may be a natural break point in the input file $inputFile");
         # result my not be undef, the assembly just has 2 or more contigs
         if ($HCov < 400) { # has cut at low coverage
            # don't try trim because there is not much to trim.
            $turnFrac = -1;
         }
      }
   }
   $log->info("Turning fraction: " . sprintf("%.4f", $turnFrac));
   # cov/frac is CPF
   return ($turnFrac, $result, $HCov/$HFrac);
}

sub lengthAdjust {
   my $len = shift;
   my $y = -0.009*$len + 2.5;  
   if ($y < 1) {
      $y = 1;
   }
   return $y; 
}

=head2 coverageToFraction 

   parameter(requiredCoverage, totalBases)
   @return (fractionOfTotal, ActualCoverage)

   When there is not enough totalBases, the 
   actual coverage will be lowered than required.
   In most cases, the requirement will be met.

=cut
sub coverageToFraction {
   my $coverage = shift;
   my $totalBases = shift;

   my $basesNeeded = int($coverage * $plasmidLength);
   if ($basesNeeded > $totalBases) { # requirement over supply
      $log->warn("coverage $coverage requires more than $totalBases bases available");
      return (1, $totalBases/$plasmidLength);
   }
   return ($basesNeeded/$totalBases, $coverage);
}

=head2 findHighestMiraFraction

parameters($inputFile, $Lower_fraction, $Higher_fraction,
   $LowerMiraResult, $HigherMiraResult

   return (bestFraction, MiraResult_of_bestFraction)

   parameter($inputFastqFile, $LowFraction, $HighFraction)
   Lower is always good, and High is always bad.

   The acend speed is controlled by M = L + (H - L)/n, now n = 3.
   The higher n, the slower the ascend speed, thus more 
   good alignment.  Maybe 4 is a good number to use.
   In Bar10 example, L starts from 0.03, H 0.616
   The next one is 0.215 too high.
   Then 0.0822

   Use STOP_CONDITION to control the termination of the
   search algorithm.

=cut
sub findHighestMiraFraction {
   my $file = shift;
   my $L = shift;
   my $H = shift;
   my $LResult = shift;
   my $HResult = shift;

   # this could be a parameter.
   if ($H - $L < $STOP_CONDITION) { # no more search
      $log->info(sprintf("%.6f", $L) . " and " . sprintf("%.6f", $H) . " close enough, no more try");
      return ($L, $LResult);
   }
   # need to do more assembly, use lower 1/3, this way
   # we do more computation at the lower end
   # We could use 1/4 if this is better.
   #my $M = $L + ($H - $L)/$ASCENT;
   #$log->debug("$M = $L + 1/3 of ( $H - $L )");
   my $M = computeMidpoint($L, $H);
   #$log->debug("Mid point $M from $L and $H");
   OUTTRAP:   my ($status, $result) = runAndEvaluateMiraFraction($file, $M);
   if ($status == 1) { # good quality
      $L = $M;
      $LResult = $result;
   }
   else {
      # M could fall in the trap, in bar16
      # 0.482 is good, 0.5 is bad
      # 0.28498 is good but 0.28945 is bad
      if (($H - $M) > 2*$STEP) {
         $log->warn("potential local trap at $M trying to get out of it with new mid point between it and $H");
         $M = ($H + $M)/2;
         goto OUTTRAP;
      }
      $H = $M;
      $HResult = $result;
   }
   return findHighestMiraFraction($file, $L, $H, $LResult, $HResult);
}

=head2 computeMidpoint

 Right now I am using linear search from the bottom.
 This is for testing purpose. I will reduce it to the
 fast search method later.

=cut
sub computeMidpoint {
   my ($L, $H) = @_;
   if ($H - $L > $STEP + $STOP_CONDITION) {
      #$log->debug("H $H L $L step $STEP (H-L) " . ($H - $L));
      return $L + $STEP; # don't step over STEP length
   }
   else {
      #return $L + ($H - $L)/$ASCENT;
      return (($ASCENT-1)*$L + $H)/$ASCENT;
   }
}

=head2 makePartialFileName

 parameter($inputFastqFileName, $fraction)

=cut
sub makePartialFileName {
   my $file = shift; # always the same file such as bar1.fastq
   my $frac = shift; # such as 0.7

   my $partfile = $file;
   my $millage = fraction2string($frac);
   if ($file =~ /\.fastq$/) {
      $partfile =~ s/\.fastq/_$millage.fastq/;
   }
   else {
      $partfile = $file . "_$millage.fastq";
   }
   #$log->debug("$millage partial file name: $partfile");
   return $partfile; # partial file
}

=head2 fraction2string

Helper function to convert fraction to 4 or less digit integer string.

=cut
sub fraction2string {
   my $frac = shift;
   #$log->debug("fraction is $frac");
   # the old method
   #return int(1000 * $frac);
   return int(sprintf("%.4f", $frac)*10000);
}

sub countFastq {
   my $file = shift;
   my $tmp = `wc -l $file`;
   if ($tmp =~ /^(\d+)/) {
      return $1/4;
   }
   else {
      die "cannot count number of sequences in file $file\n";
   }
}

=head2 evaluateFastq

return (mean, count, total) bases in longth

=cut
sub evaluateFastq {
   my $fastqFile = shift;
   my $statFile = getFileStem($fastqFile) . ".stat";
   my ($count, $mean, $total);
   if (!-f $statFile || -z $statFile) {
      system("$binaryDir/fastqstat $fastqFile $statFile");
   }
   # now we should have fastq file for sure
   open IN, "<$statFile" || die "Failed to open $statFile for reading!\n";

   # read the simple one line stat file
   my $header = <IN>;
   chomp $header;
   my @hr = split /\t/, $header;
   my $line = <IN>;
   chomp $line;
   if ($hr[0] eq 'mean' && $hr[2] eq 'count') {
      my @row = split /\t/, $line;
      $count = $row[2];
      $mean = $row[0];
   }
   else { # other format
      $line = <IN>;
      while ($line && $line !~ /after/) {
         $line = <IN>;
      }
      chomp $line;
      my @row = split /\t/, $line;
      $count = $row[3];
      $mean = $row[1];
   }
   $total = sprintf("%i", $mean*$count);

   close IN;
   return ($mean, $count, $total);
}


=head2 trimsequences

   Trim the 3'-low-quality part from the reads using the trimfastq program.

   parameter(inputfile)
   return reference to an array of trimmed files that are determined by the
      trim parameter array.

=cut
sub trimsequences {
    my $infile = shift;

    my @trimmed;
    my $doneFile = $infile . "_trim.done";
    if (-f $doneFile) {
        $log->debug("$infile trim has been done before, not repeating");
        my $stem = getFileStem($infile);
        my $filePattern = $stem . "_trim[3-6][0-9][0-9].fastq";
        # if file pattern failed we can use the done file to store the results.
        @trimmed=`ls $filePattern`;
        if (!scalar(@trimmed)) {
           die "there is only done file but no trim result!\n";
        }
        chomp @trimmed;
        # order the file in the order 420,419 520,519, 620, 619
        for (my $i = 0; $i < @trimmed-1; $i += 2) {
           if ($trimmed[$i] =~ /[4,5,6]19/ && $trimmed[$i+1] =~ /[4,5,6]20/) {
              my $x = $trimmed[$i];
              $trimmed[$i] = $trimmed[$i+1];
              $trimmed[$i+1] = $x;
           }
           else {
              $log->warn("trim file pattern strange " . $trimmed[$i] . " and "
                 . $trimmed[$i+1]);
           }
        }
        return \@trimmed;
    }
    my $outfile;

    foreach my $param (@trimParameters) {
        my $info = `trimfastq -i $infile -w $param->[0] -c $param->[1]`;
        if ($info =~ /written to (.+)$/m) {
            $outfile = $1;
        }
        else {
            die "$info has no output file name\n";
        }
        #print "result written to $outfile\n";
        push @trimmed, $outfile;
    }
    system("touch $doneFile");
    return \@trimmed;
}

=head2 linearizeInput
   
   linearize the input reads in the file.
   parameter(inputfile)
   return the linearized file name.

=cut
sub linearizeInput {
   my $file = shift;
   my $cutfile = getFileStem($file) . "cut.fastq";
   my $doneFile = $cutfile . ".done";
   if (-f $doneFile) {
      $log->info("cutfile $cutfile done before, saw donefile: $doneFile, not repeating");
      return $cutfile;
   }
   else {
      $log->info("linearizing $file ...");
      system("$binaryDir/linearize $file $cutfile");
      if ($?>8) {
         $log->fatal("Failed to linearize $file while running $binaryDir/linearize $file $cutfile");
         #$log->fatal("Failed to run linearization on $file");
         die "Failed to run linearization on $file";
      }
      system("touch $doneFile");
      $log->info("linearized result file: $cutfile");
   }
   return $cutfile;
}

sub trimAndLinearize {
   my $infile = shift;
   my $trimmed = trimsequences($infile);
   my @linear;
   foreach my $t (@$trimmed) {
      push @linear, linearizeInput($t);
   }
   return \@linear;
}


=head2 runMira

   Run one mira job for one config file.

parameter(input_file)

If fail then this method will die. Mira alway run successfully regardless of
input. Unless we see other behavors. This version 4.0.  The behavior of
mira has not been stable in the past. The author is kind not sticking to
a particular interface. This cause some pain and rewrite for me.

return 1 run success but result needs to be evaluated. 
       0 for bad, cannot be evaluated, means too much input.

=cut
sub runMira {
   my $inputFile = shift;

   my $conf = fromInputToConfigFile($inputFile);
   $log->info("run mira with config: $conf ...");
   my $doneFile = $conf . ".done";
   my $logfile = $conf . ".log";
   my $status = readDoneFile($doneFile);
   if ($status != -1) {
      $log->debug("mira job $conf done before with status $status, not repeating");
      return $status;
   }

   $log->debug("A fresh run of mira (may take a long time) with $conf");
   writeMiraConfigFile($inputFile);
   # the mira binary should be controlled
   #system("mira $conf > $logfile");
   system("$miraBinary $conf > $logfile");
   if ($?>>8) {
      my $reason = getMiraErrorFromLogFile($logfile, $conf);
      if ($reason) {
         if ($reason ==  9) {
            die "Mira died from some unknown reason!\n";
         }
         else {
            die "Mira died because of unrecoverable problem\n";
         }
      }
      writeDoneFile($doneFile, 0);
      return 0; # failure too much input
   }
   writeDoneFile($doneFile, 1);
   print STDERR "mira log file: $logfile has more detailed information.\n";
   return 1;
}

=head getMiraErrorFromLogFile

 return 1 for duplicated sequence names
        2 for disk full
        0 for too much input, just could not handle
        9 for other errors that requires manual examination of the
          log file.

=cut
sub getMiraErrorFromLogFile {
   my $logfile = shift;
   my $conf = shift;

   open LO, "<$logfile" or die "Could not open logfile: $logfile!\n";
   my $line = <LO>;
   while ($line) {
      if ($line =~ /Could not write anymore to disk/) {
         $log->fatal("Disc full");
         return 2;
      }
      elsif ($line =~ /maximum ratio has been reached/) {
         $log->warn("Failed run because too much input for mira, try less?");
         return 0;
      }
      elsif ($line =~ /Some read names were found more than once/) {
         $log->fatal("duplicated sequence name, fix it first before running it again.");
         return 1;
      }
      $line = <LO>;
   }
   $log->fatal("Failed to run with $conf. Please check $logfile");

   return 9;
}

=head2 evaluateMira

return an array the first element: 1 for good 0 for bad
    the second element is the MiraResult
    if there is no mira result, then it is undef.

    Even the first element could be zero, the second element may contain
    mira result that are bad.

=cut
sub evaluateMira {
   my $inputFile = shift;

   $log->debug("evaluating mira with input: $inputFile ...");
   my $stem = getFileStem($inputFile); # project name
   my $miraResult = new MiraResult($stem);
   return ($miraResult->evaluateResult, $miraResult);
}

=head2 runAndEvaluateMira

   return an array of two elements. 
      The first element:1 for good 0 for bad, 2 for two large contigs
      The second element is the MiraResult.

=cut
sub runAndEvaluateMira {
   my $inputFile = shift;

   # runMira should return 1 for good, 0 for bad, 2 for too many
   # use the same status as evaluate Mira.
   my $status = runMira($inputFile);
   if ($status == 0) { # too much input result 0 status
      # and cannot be evaluated because no output produced.
      return (0, undef);
   }
   return evaluateMira($inputFile);
}

=head2 runAndEvaluateMiraFraction

parameters($inputFileName, $fraction)

return an array of (the status of the run, MiraResult)
   0 status for too much input or zero large contigs or too many contigs

make runAndEvaluateMira a special case of runAndEvaluateMiraFraction

=cut
sub runAndEvaluateMiraFraction {
   my $file = shift;
   my $frac = shift;

   $log->info("run and evaluate $file with fraction " . sprintf("%.5f", $frac) . " ...");
   if ($frac == 1) {
      return runAndEvaluateMira($file);
   }
   else {
      my $partialFile = makePartialFileName($file, $frac);
      $log->debug("Trying to make a fraction fastq input file: $partialFile ...");
      my $cmd = "$binaryDir/cpfasqpart -p $frac -o $partialFile $file";
      system($cmd);
      if ($?>>8) {
         $log->fatal("Failed to run $cmd!");
         die "Failed to run $cmd\n";
      }
      return  runAndEvaluateMira($partialFile);
   }
}

sub writeDoneFile {
   my $fileName = shift;
   my $status = shift;

   open OU, ">$fileName" or die "Cannot open $fileName for writting\n";
   print OU "status=$status\n";
   close OU;
}

=head1 readDoneFile

 return -1 if no such file. otherwise return the status of the job done

=cut
sub readDoneFile {
   my $fileName = shift;
   if (open(DONE, "<$fileName")) {
      my $status = <DONE>;
      chomp $status;
      if ($status =~ /status=(\d+)/) {
         close DONE;
         return $1;
      }
      else {
         close DONE;
         die "unexpected status line $status\n";
      }
   }
   else {
      if (!-f $fileName) {
         return -1;
      }
      else {
         warn "Cannot open done file $fileName even it is there!\n";
         return -1;
      }
   }
}

=head2 writeMiraConfigFile

   parameter(inputfile)
   Use one file per input. Mira can take multiple files,
   but we are not using that for better automation control.

   return config file name. The config file is made of the stem
          of the input file with the conf suffix.

   The plasmids are in the range of 9 kb.

=cut
sub writeMiraConfigFile {
    my $fileName = shift;
    my $stem = getFileStem($fileName);
    my $configFile = fromInputToConfigFile($fileName);
    my $largeContigForStat = int($plasmidLength/3);
    # rewrite it anyway, in case I need to modify the config
    #if (-f $configFile) {
    #    $log->debug("configFile exists not writting $configFile");
    #    return $configFile;
    #}

    my $content=<<ENDS;
project = $stem
job = genome,denovo,accurate
readgroup = Shotgun
data = $fileName
technology = iontor
parameters = COMMON_SETTINGS -GENERAL:number_of_threads=4 -AS:automatic_repeat_detection=off -AS:uniform_read_distribution=off -NW:check_nfs=no -DI:tmp_redirected_to=$tmpdir -MI:large_contig_size=450 -MI:large_contig_size_for_stats=$largeContigForStat IONTOR_SETTINGS -AS:mrl=17 -AS:coverage_threshold=6.5
ENDS
    # -MI:sonfs=no  seem to have be disabled or moved to -NW, 
    # The Mira documentation is still there.
#parameters = COMMON_SETTINGS -AS:automatic_repeat_detection=off -AS:uniform_read_distribution=off -MI:large_contig_size=450 -MI:large_contig_size_for_stats=3000 IONTOR_SETTINGS -AS:mrl=17 -AS:coverage_threshold=3.5 -AL:min_overlap=13 
    # if min_overlap=13, then linear cut will not be effective,
    # removing this 
    open OU, ">$configFile" or die $!;
    print OU $content;
    close OU;
    return $configFile;
}

=head2 getFileStem

I am using my own function over the perl standard
because of the simple nature of this method. And I may
have better control.

File::Basename has all the features
($name, $dir, $ext) = fileparse($path, '\..*');

=cut
sub getFileStem {
   my $fileName = shift;
   my $stem;
   if ($fileName =~ /(.+)\.[a-z]+$/) {
      $stem = $1;
   }
   else { $stem = $fileName; }
   return $stem;
}

sub fromConfigToInputFile {
   my $config = shift;
   my $input = getFileStem($config)  . ".fastq";
   return $input;
}

sub fromInputToConfigFile {
   my $input = shift;
   my $config = getFileStem($input) . ".conf";
   return $config;

}

sub getLogFileName {
   #if (!$ENV{HOME}) { die "no home\n"; }
   if (!$inputseq) { die "no inputseq\n"; }
   # this can be a problem
   #return $ENV{HOME} . "/log/" . $inputseq . ".log";
   # use current directory for log
   return $inputseq . ".mirarun.log";
}

sub updateStatusDone {
   my $datasource = 'dbi:Oracle:bioinfo';
   my $dbh = DBI->connect($datasource, 'miramill', 'miramill', { RaiseError => 1, AutoCommit => 0 })
      or die "Connection to database $datasource failed: ", $DBI::errstr, "\n";
   #$dbh->{LongReadLen}=60372;
   my $projname = $inputseq;
   $projname =~ s/\.fastq$//;
   my $sth = $dbh->prepare("update request set projstatus=(select id from project_status where description='assembled') where name='$projname'");
   eval {
      $sth->execute;
   };
   if ($@) {
      $log->fatal("Failed to upload project status in database: " . $dbh->errstr);
      $dbh->rollback or die $dbh->errstr;
      return 0;
   }
   else {
      $dbh->commit or die $dbh->errstr;
      $log->info("project status updated to 'compared to ref'");
      return 1;
   }
}

__END__

=head1 NAME

runmira - simple pipeline to assembly plasmid

=head1 SYNOPSIS

 runmira bar1.fastq

 when the above is successful
 runmira -c in the current run directory to compare to the reference.
 The above run depends on 1. successful mira run, and 2. has refseq

=head1 DESCRIPTION

Assembly plamid with mira from ion torrent sequence reads.
One input file is used.  All the work will be done in one directory.
Dependign on the number of input sequences, this prodedure will do trimming and linearization.
Trimming will get rid of noisy information only if the input sequence is overwhelming.
Linearization will be done if there is no natural breaks in the assembly.

=head2 Options

 -l plasmid length. This can also be given in the configuration file.
 -H HIGH_COVERAGE default 1500
 -L LOW_COVERAGE default 100
 --flip-consensus   flip the consensus of assembly when comparing assembly with reference.
 --flip-reference   flip the referece when comparing to assembly
 -A ASCENT  default 3. Orginary user should not use this option.
    This parameter controls the trajectory of the midpoint calculation.
    The larger the slower the climing.  2 would be the fastest.
 --clean clean the running directory to start a new run.
 -c  flag to only check the quality of the latest run. You must have
     successfully run a mira before using this command, and you also
     have to have the reference sequence. 
 --help print help message.

=head2 PROBLEMS

When the plasmids are too short, the starting input sequence should be lowered.

=head1 Requirements

To run this program, you need to follow the following rules and input files:

=head2 Input file

input file must be named: bar.fastq

This program will only take fastq as input, and only one file.

=head2 config file

runmira.cfg is the default configuration file for running mira.
You should have one in your running directory.  If not present all 
of the parameters can also be given through the command line.
Here is an example file:

 input=bar13.fastq
 length=6372
 low=100
 high=1400
 

=head2 linearize.cfg

This file contains the linearization cut-site sequence and other parameters.
You must have this file.

example file content

   head=CTACCAGCGGTGGTTTGTTTGCCGGATCAAGAGCTA
   tail=CGCGTAATCTGCTGCTTGCAAACAAAAAAACCACCG
   identity-cutoff=0.84
   alnlen-cutoff=17
   iden-discard=0.92
   gap-discard=2
   shortest=16

