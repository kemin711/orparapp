package MiraResult;

use strict;
use Log::Log4perl;
use FastaReader;

=head1 NAME

MiraResult - a quick and dirty representation of mira assembly result

=head1 SYNOPSIS

my $mr = new MiraResult($fileStem)
   The basename of the file (or stem) such as the bar1.fastq file
   you need to use bar1 as the parameter for this function.

 $mr->evaluateResult;

=head2 new

Constructor($fileStem)

   FileStem is also the project name.
   Mira use this to build up path and file names.

   This object is only meant to be used in the current directory.
   It uses the current directory as working directory.

=cut
sub new {
   my $invocant = shift;
   my $stem = shift;

   my $class = ref $invocant || $invocant;
   my $inputFile = $stem . ".fastq";
   my $numseq = `wc -l $inputFile`;
   if ($numseq =~ /^(\d+)/) {
      $numseq = $1/4;
   }
   else {
      die "Cannot count number of sequences in file $inputFile!\n";
   }
   my $self = {proj => $stem, numreads => $numseq };
  
   bless $self, $class;
   $self->readInfo;
   return $self;
}

=head2 getInputFile

Right now fixed to one input file per project.
These are the raw reads.

=cut
sub getInputFile {
   my $self = shift;
   return $self->getProjectName . ".fastq";
}

sub getProjectName {
   my $self = shift;
   return $self->{proj};
}

=head2 isBadResult

This method will look for the result file without fully parsing the outputs 
from mira run.

Must use this method in the directory where mira was run.

=cut
sub isBadResult {
   my $self = shift;
   my $dirname = $self->getAssemblyDirectory;
   my $results = `find $dirname -name *_info_assembly.txt | xargs grep 'Number of contigs:'`;
   my @res = split /\n/, $results;
   if (@res != 2) {
      die "expecting only two elements!\n";
   }
   my @numcontig;
   if ($res[0] =~ /(\d+)/) {
      $numcontig[0] = $1;
   }
   else {
      die "failed to extract number of largest contigs!\n";
   }
   if ($res[1] =~ /(\d+)/) {
      $numcontig[1] = $1;
   }
   else {
      die "failed to extract number of contigs!\n";
   }
   print "number of largest contigs: $numcontig[0]  all contigs: $numcontig[1]\n";
   if ($numcontig[0] == 0) {
      return 1;
   }
   else {
      return 0;
   }
}

sub getAssemblyDirectory {
   my $self = shift;
   return $self->getProjectName . "_assembly";
}

=head2 getResultPath

return the assembly result file path

=cut
sub getResultPath {
   my $self = shift;
   my $path = $self->getAssemblyDirectory . "/" . $self->getProjectName . "_d_results";
   return $path;
}

=head2 getInfoPath

return the assembly information unix file path

=cut
sub getInfoPath {
   my $self = shift;
   my $path = $self->getAssemblyDirectory . "/" . $self->getProjectName . "_d_info";
   return $path;
}

=head1 readAssemblyInfo

   reads_assembled => number
   avg_coverage => float
   longest_contig => integer
   max_coverage => integer
   quality => integer

   We only process the large contigs.

=cut
sub readAssembleInfo {
   my $self = shift;
   my $path = $self->getInfoPath . "/" . $self->getProjectName . "_info_assembly.txt";
   open IN, "<$path" or die "Cannot open $path!\n";
   my $line = <IN>;
   my %result = ();
   # the result is divided into two sections: large and all
   # large is defined as > 5000 nt long contigs.
   # This parameter is controlled by mira =MI:large_contig_size
   # stop parsing after the largest contig
   # We only collect the large contigs
   my $logger = Log::Log4perl->get_logger("MiraResult");
   while ($line && $line !~ /^All contigs:/) {
      if ($line =~ /^Num\. reads assembled:\s+(\d+)/) {
         $result{reads_assembled} = $1;
         #$logger->debug("assembled reads $1");
      }
      elsif ($line =~ /Avg\. total coverage:\s+(\d+\.\d+)/) {
         $result{avg_coverage} = $1;
         #$logger->debug("Average total coverage: $1");
      }
      elsif ($line =~ /Largest contig:\s+(\d+)/) {
         $result{longest_contig} = $1;
         $logger->debug("Longest contig: $1");
      }
      elsif ($line =~ /Max coverage \(total\):\s+(\d+)/) {
         $result{max_coverage} = $1;
         #$logger->debug("Maximum coverage: $1");
      }
      elsif ($line =~ /Average consensus quality:\s+(\d+)/) {
         $result{quality} = $1;
         #$logger->debug("Average consensus quality: $1");
      }
      elsif ($line =~ /Number of contigs:\s+(\d+)/) {
         $result{numcontig} = $1;
         $logger->debug("Number of large contigs: $1");
      }
      $line = <IN>;
   }
   $self->{info} = \%result;
   close IN;
   return \%result;
}

sub getAssembledReads {
   my $self = shift;
   return $self->{info}{reads_assembled};
}

sub getCoverage {
   my $self = shift;
   return $self->{info}{avg_coverage};
}

=head2 getLongestContig
   
 return the length of longest contig

=cut
sub getLongestContig {
   my $self = shift;
   return $self->{info}{longest_contig};
}

sub getMaxCoverage {
   my $self = shift;
   return $self->{info}{max_coverage};
}

=head2 getQuality

 return the average consensus quality

=cut
sub getQuality {
   my $self = shift;
   return $self->{info}{quality};
}

sub getNumberOfContigs {
   my $self = shift;
   return $self->{info}{numcontig};
}

=head readAssembleStat

   The original rows from the mira assembler
   The following is the header for the bar3cut_info_contigstats.txt
   # name length   av.qual  #-reads  mx.cov.  av.cov   GC% CnIUPAC  CnFunny  CnN   CnX   CnGap CnNoCov
   To me the most important fields are:
      0: name 1: av.qual (average quality) 2:#-reads 4:av.cov

   It is not worth to make an object for this at this point.
   May be in the future.

=cut
sub readAssembleStat {
   my $self = shift;
   my $path = $self->getInfoPath . "/" . $self->getProjectName . "_info_contigstats.txt";
   open IN, "<$path" or die "Cannot open $path\n";
   my $line = <IN>;
   chomp $line;
   my @header = split /\t/, $line;
   $self->{contigstat_header} = \@header;
   my @table = ();
   while (<IN>) {
      chomp;
      my @row = split /\t/;
      @row = map { s/\s+$//; $_; } @row;
      push @table, [ @row ];
   }
   close IN;
   $self->{contigstat} = \@table;
   return \@table;
}

sub getAssembleStat {
   my $self = shift;
   return $self->{contigstat};
}

=head2 getFastaResultFile

return the file path (relative to the current directory)
for the unpadded.fasta file. It is usually 
   /*_assembly/*_d_results/*_out.unpadded.fasta

=cut
sub getFastaResultFile {
   my $self = shift;
   my $fasfile = $self->getResultPath . "/" . $self->getProjectName . "_out.unpadded.fasta";
   return $fasfile;
}

=head2 getQualityFile

 return the quality file name.  
 example: bar9_trim420cut_7493_out.unpadded.fasta.qual

=cut
sub getQualityFile {
   my $self = shift;
   return $self->getFastaResultFile . ".qual";
}

=head2 getFastaSequenceById

   parameter($seqid)
return a Fastaseq object from the unpadded fasta assembly result.

=cut
sub getFastaSequenceById {
   my $self = shift;
   my $seqid = shift;
   my $reader = new FastaReader($self->getFastaResultFile);
   return $reader->getByName($seqid);
}

=head2 getQualityById

 parameter($seqid)
 return the quality as a reference to an array

=cut
sub getQualityById {
   my $self = shift;
   my $seqid = shift;
   my $qfile = $self->getQualityFile;
   open QUA, "<$qfile" or die "Failed to open quality file: $qfile $!\n";
   my $line = <QUA>;
   while ($line && $line !~ /$seqid/) {
      $line = <QUA>;
   }
   if (!$line) {
      die "failed to find $seqid from $qfile\n";
   }
   $line = <QUA>;
   my @quality = ();
   while ($line && $line !~ /^>/) {
      chomp $line;
      my @arr = split /\s+/, $line;
      push @quality, @arr;
      $line = <QUA>;
   }
   return \@quality;
}

=head2 readInfo

 this is the initialization method.
 I will put this into the new method.

=cut
sub readInfo {
   my $self = shift;
   if (!$self->{initialized}) {
      my $inputFile = $self->getInputFile;
      my $numseq = `wc -l $inputFile`;
      if ($numseq =~ /^(\d+)/) {
         $numseq = $1/4;
      }
      else {
         die "Cannot cound number of sequences in file $inputFile!\n";
      }
      $self->{numreads} = $numseq;
      $self->readAssembleInfo;
      $self->readAssembleStat;
      # use state variable, not good design
      $self->{initialized} = 1;
   }
}

sub getNumberOfReads {
   my $self = shift;
   return $self->{numreads};
}

sub getAssembledReadFraction {
   my $self = shift;
   if ($self->getNumberOfReads == 0) {
      die $self->getProjectName, " has zero number of reads!, may be crashed run.\n";
   }
   return $self->getAssembledReads / $self->getNumberOfReads;
}

=head2 getGoodContigs

Actually is getting the names of the good contigs
not the actual object.  This is just a quick and dirty work.

=cut
sub getGoodContigs {
   my $self = shift;
   my $logger=Log::Log4perl->get_logger("MiraResult");
   
   if (exists $self->{goodContigs}) {
      return $self->{goodContigs};
   }

   my @goodContigs;
   my $i=0;
   foreach my $s (@{$self->getAssembleStat}) {
      if ($s->[5] > $self->getCoverage*0.8) {
         ++$i;
         push @goodContigs, $s->[0];
      }
   }
   $self->{goodContigs} = \@goodContigs;
   return \@goodContigs;
}

sub getNumberOfGoodContigs {
   my $self = shift;
   return scalar(@{$self->getGoodContigs});
}

=head2 getLargestContigSequence

return the Fastaseq object with for the largest contig.

=cut
sub getLargestContigSequence {
   my $self=shift;
   my $logger=Log::Log4perl->get_logger("MiraResult");
   my $goodid = $self->getGoodContigs;
   if (!$goodid) {
      $logger->warn("There is no good contig in this assembly!");
      return undef;
   }
   return $self->getFastaSequenceById($goodid->[0]);
}

=head2 evaluateResult

return 1 for good (only one contig), 2 for having two contigs, 0 for bad.

Good means there is only and exactly one long contig

=cut
sub evaluateResult {
   my $self = shift;
   if (!$self->{initialized}) {
      $self->readInfo;
   }
   my $logger=Log::Log4perl->get_logger("MiraResult");

   my $good = 0; # default is bad
   #$logger->debug("Fraction of reads assembled: " .  sprintf("%.3f", $self->getAssembledReadFraction));
   $logger->debug(" average coverage of large contigs: " . $self->getCoverage);
   #$logger->debug("good contigs:");
   if ($self->getAssembledReadFraction < 0.3 || $self->getCoverage < 5) {
      $logger->warn("Something is wrong with the assembly, and evaluation is not meaningful!");
      $good = 0;
   }
   my $goodcontig = $self->getGoodContigs;

   if (@$goodcontig > 18) {
      $logger->warn("There are " . scalar(@$goodcontig) . " good contigs, too many!");
      $good = 0;
   }
   elsif (@$goodcontig == 1) {
      #$logger->debug("there is only one large contig, result good");
      return 1;
   }
   elsif (@$goodcontig == 2) {
      #$logger->debug("There are two good contigs: " . join(" | ", @$goodcontig));
      return 2;
   }

   if ($self->getLongestContig == 0) {
      $logger->warn("there is no longest contig, bad");
      $good = 0;
   }
   $logger->info("evaluation result: $good");
   return $good;
}

=head2 showBasicInformation

print the length of longest contig, consensus quality, average coverage

=cut
sub showBasicInformation {
   my $self = shift;
   my $fh = shift;
   if (!$fh) { $fh = \*STDOUT; }
   print $fh "length largest contig: ", $self->getLongestContig, ", ",
         "consensus quality: ", $self->getQuality,  ", ",
         "coverage: ", $self->getCoverage, "\n";
}

=head2 showBreakPoint

If the assembly is done without linearization then there is a natural
break point in the input sequence.

=cut
sub showBreakPoint {
   my $self = shift;
   my $fh = shift;
   if (!$fh) { $fh = \*STDOUT; }

   my $fastaSeq = $self->getLargestContigSequence;
   my $head = $fastaSeq->getFirst(40);
   #print "first\n";
   #$head->show;
   my $tail = $fastaSeq->getLast(40);
   #print "last\n";
   #$tail->show;

   print $fh "break point:\n",
      $tail->sequence, " | ",  $head->sequence, "\n";

   my $headrc = $head->revcompCopy();
   my $tailrc = $tail->revcompCopy();
   print $fh "break viwed from the reverse complement\n",
      $headrc->sequence, " | ", $tailrc->sequence, "\n";
}

=head2 getEndSequences

   return reference to an array (head, tail, head-rc, tail-rc)

=cut
sub getEndSequences {
   my $self = shift;
   my $length = shift;

   if (!$length) { $length = 40; }


   my $fastaSeq = $self->getLargestContigSequence;
   my $head = $fastaSeq->getFirst($length);
   my $tail = $fastaSeq->getLast($length);

   my $headrc = $head->revcompCopy();
   my $tailrc = $tail->revcompCopy();
   return [$head->sequence, $tail->sequence, $headrc->sequence, $tailrc->sequence];
}

1;
