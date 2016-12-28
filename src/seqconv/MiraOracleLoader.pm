package MiraOracleLoader;

## Loader of mira results ###

# to generate a table of request to be loaded into the database
# This program has to be run from the directory contain the mira runs

use strict;
use warnings;
use FastaReader;
use Fastaseq;
use DBI;
use Log::Log4perl;

my $log = Log::Log4perl->get_logger("MiraOracleLoader");

=head1 NAME

 MiraOracleLoader - loading mira run results with reference sequence to oracle database

=head1 SYNOPSIS

   my $loader = MiraOracleLoader->new('bar10');
   $loader->setRunDirectory('/net/kraken/ng10/kzassembly');
   $loader->upload;

=head2 new

  Use this module in the top directory
  This hard-coded on the production database.
  All work will be on production server.

  The underlying connection will not do autocommit.
  This may save us from connection problems.

  The database is the production version on Bioinfo.

=cut
sub new {
   my $invocant = shift;
   my $projname = shift;
   if (!$projname) {
      die "you must provide a project name\n";
   }

   my $class = ref $invocant || $invocant;
   $log->debug("connecting to the miramill\@bioinfo database");
   my $dbh = DBI->connect('dbi:Oracle:bioinfo', 'miramill', 'miramill', 
         { RaiseError => 1, AutoCommit => 0 })
         or die "Connection to bioinfo failed: ", $DBI::errstr, "\n";
   $log->debug("setting up FastaReader without input file");
   my $seqreader = new FastaReader();

   my $hashref = { 
      dbh => $dbh,
      #rundir => '/d1/work/assembly',
      rundir => '/data/assembly',
      projname => $projname, 
      #requester => 'Xin Huang',
      #description => 'Mira automatic assembly',
      outfile => 'miraproj.tab',
      resultdir => 'results',
      qcfile => 'QC.tab',
      resultfile => 'bestResult.tab',
      seqreader => $seqreader
   };
   bless $hashref, $class;
   return $hashref;
}

sub setProjectName {
   my $self = shift;
   $self->{projname} = $_[0];
}

=head2 getProjectName

 return the project name

=cut
sub getProjectName {
   my $self = shift;
   return $self->{projname};
}

sub getQCFileName {
   my $self = shift;
   return $self->{qcfile};
}

sub getResultDirectory {
   my $self = shift;
   return $self->{resultdir};
}

sub getResultFileName {
   my $self->shift;
   return $self->{resultfile};
}

sub getQC {
   my $self = shift;
   return $self->{qc};
}


=head2 fetchProjectInfo

 fetch the project information from the database based on the project name.
 project name (name column for table request) is a unique key.

=cut
sub fetchProjectInfo {
   my $self = shift;

   my $projname = shift;
   if (!$projname) {
      if (!$self->{projname}) {
         die "you either provide a project name here or at object construction time\n";
      }
      $projname = $self->{projname};
   }

   $log->debug("fetching data for $projname ...");

   $self->{dbh}->{LongReadLen}=500000;
   my $sth = $self->{dbh}->prepare("select id, description, head_len, tail_len, requester, refseq from request where name='$projname'");
   $sth->execute;
   my @row = $sth->fetchrow_array;
   if ($sth->rows == 0) {
      $log->debug(join(' | ', @row));
      die "you have not entered the project info through the web page, do it before you run upload\n";
   }
   #$log->debug("content for fetch: ", join(" | ", @row));
   my %info = ( id => $row[0], name => $projname, description => $row[1],
      head_len => $row[2], tail_len => $row[3], requester => $row[4],
      refseq => $row[5] );
   $self->{projinfo} = \%info;
   return \%info;
}

sub showProjectInfo {
   my $self = shift;
   $log->info("project: ", $self->getProjectName);
   my $info = $self->{projinfo};
   foreach my $k (keys %$info) {
      print "$k => $info->{$k}\n";
   }
}

sub getRequestId {
   my $self = shift;
   if (!$self->{projinfo}) {
      die "you have not downloaded project informaton for ", $self->getProjectName, "\n";
   }
   return $self->{projinfo}->{id};
}

# Disable modification of database data
#sub setDescription {
#   my $self = shift;
#   $self->{description} = $_[0];
#}

=head2 getDescription
 
 obtain the description field for the project

=cut
sub getDescription {
   my $self = shift;
   if (!exists $self->{projinfo}) {
      die "project info not loaded yet!\n";
   }
   return $self->{projinfo}{description};
}

# disable modification of database data
#sub setRequester {
#   my $self = shift;
#   $self->{requester} = $_[0];
#}

sub getRequester {
   my $self = shift;
   if (!exists $self->{projinfo}) {
      die "no project loaded!\n";
   }
   return $self->{projinfo}{requester};
}

sub setProjectDirectory {
   my $self = shift;
   $self->{workdir} = $_[0];
}

=head2 setRunDirectory
  
  set up the directory where mira run was done.
  This is sort of the home directory for mira run.

=cut
sub setRunDirectory {
   my $self=shift;
   if (! -d $_[0]) {
      die "directory $_[0] does not exists\n";
   }
   $self->{rundir} = $_[0];
}

sub getRunDirectory {
   my $self=shift;
   return $self->{rundir};
}


=head2 getRefseqFile

 The default refseq file name is named in the work directory
 projname/projnameknown.fas
 This is for the convenice of constructing the pipeline.
 We are following certain conventions.

 return a full path to the referencce sequence fasta file.

=cut
sub getRefseqFile {
   my $self = shift;
   return $self->getRunDirectory . $self->getProjectName . "/" . $self->getProjectName . 'known.fas';
}

=head2 getWorkDirectory

 The work directory is where the mira is run. All result files are stored here.
 It is typically MIRA_ROOT/projname

=cut
sub getWorkDirectory {
   my $self = shift;
   return $self->getRunDirectory . '/' . $self->{projname};
}

sub getResultfile {
   my $self = shift;
   return $self->{resultfile};
}

sub getSequenceReader {
   my $self = shift;
   return $self->{seqreader};
}

=head2 upload

 This method is intended to be used by users. 
 It upload the result into the database.

 First upload result, then QC info.

 returns true if success false otherwise

=cut
sub upload {
   my $self = shift;
   $self->fetchProjectInfo;
   if ($self->uploadResult) {
      if ($self->uploadQC) {
         return 1;
      }
      else {
         $log->error("faile to upload QC");
         return 0;
      }
   }
   else {
      $log->error("Failed to upload result");
      return 0;
   }
}

=head2 markDone

 Run this method if project is successfully done.
 This method will update the project status.

=cut
sub markDone {
   my $self = shift;

   my $projname = $self->getProjectName;
   my $dbh = $self->{dbh};
   my $sth = $dbh->prepare("update request set projstatus=(select id from project_status where description='compared to ref') where name='$projname'");

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

####################################################
####################### qc data upload #############
# all work is done in the project subdirectory
sub uploadQC {
   my $self = shift;
   $log->info("uploadingQC results to database");

   my $workdir = $self->getWorkDirectory;
   my $pwd = `pwd`;
   chomp $pwd;
   $log->debug("working directory is $workdir. And we are in $pwd");
   if ($pwd ne $workdir) {
      chdir $workdir;
   }
   eval {
      $self->readQC;
      $self->transformQC;
      if (!$self->uploadQCToDatabase) {
         return 0;
      }
      if (!$self->integrateQC) {
         return 0;
      }
   };
   if ($@) {
      $log->fatal("Failed uploadQC: $@");
      return 0;
   }
   else {
      chdir "-";
      $log->info("Uploaded QC to database done");
      return 1;
   }
}

sub integrateQC {
   my $self = shift;
   $log->info("integrateQC");

   my $dbh = $self->{dbh};
   my $sqlstr=<<ENDQ;
insert into compare_ref
select r.id,
   c.ends, c.direction_ref, c.direction_ass, c.alnlength, 
   c.numgap_ref, c.numgap_ass, c.gaplen_ref, c.gaplen_ass, 
   c.identity, c.begin_ref, c.end_ref, c.begin_ass, c.end_ass,
   c.alntext
from raw_compare_ref c join result r on c.input = r.input and c.fraction = r.fraction
ENDQ
   eval {
      my $sth = $dbh->prepare($sqlstr);
      $sth->execute or die $dbh->errstr;
   };
   if ($@) {
      $log->fatal("Failed to integrate QC result");
      $dbh->rollback or die $dbh->errstr;
      return 0;
   }
   else {
      $dbh->commit or die $dbh->errstr;;
      return 1;
   }
}

=head2 uploadQCToDatabase

 return true for success false for failure

=cut
sub uploadQCToDatabase {
   my $self = shift;
   $log->info("uploadQCToDatabase");

   my $qc = $self->getQC;
   my $dbh = $self->{dbh};
   my $request_id = $self->getRequestId;
   # before loading into staging table, the staging table needs to be cleared.
   my $sth = $dbh->prepare("truncate table raw_compare_ref");
   $sth->execute or die $dbh->errstr;
   $dbh->commit or die $dbh->errstr;

   # load into staging table
   $sth = $dbh->prepare("insert into raw_compare_ref (input, fraction, ends, direction_ref, direction_ass, alnlength, numgap_ref, numgap_ass, gaplen_ref, gaplen_ass, identity, begin_ref, end_ref, begin_ass, end_ass, alntext) values (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)");

   my $h = $qc->[0];
   $log->debug("result header:\n" . join("\t", @$h));
   eval {
      for (my $i=1; $i < @$qc; ++$i) {
         $sth->execute(@{$qc->[$i]});
      }
   };
   if ($@) {
      $log->fatal("Failed to upload QC data into database: " . $dbh->errstr);
      $dbh->rollback or die $dbh->errstr;
      return 0;
   }
   else {
      $dbh->commit or die $dbh->errstr;
      $log->info("QC loaded to database raw_compare_ref");
      return 1;
   }
}

sub transformQC {
   my $self = shift;
   $log->debug("Starting transformQC");

   my $resultDir = $self->getResultDirectory;
   my $result = $self->{qc}; # qcresults

   my $h = $result->[0];
   $log->debug("headers:\n" . join("\t", @$h));
   for (my $i=1; $i < @$result; ++$i) {
      my $fileName;
      if ($result->[$i][1] == 1) {
         $fileName = $result->[$i][0] . "_c*.aln";
      }
      else {
         $fileName = $result->[$i][0] . '_' . fraction2string($result->[$i][1]) . "_*.aln";
      }
      $log->debug("Looking for *aln files $fileName in directory $resultDir");
      my @files = glob("$resultDir/$fileName");
      my $correctFile;
      if (@files == 1) {
         $correctFile = $files[0];
         # should have only one file!
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
            $log->warn("Conflicting results: " . join(' : ', @{$result->[$i]}));
            $log->warn("two aligns in different direction aligns to the reference! picking one in random: " . join(' | ', @files));
            $correctFile = $files[0];
            #die "Found two align file one is bad and the other is good\n",
            #   join(' | ', @files), "\n";
         }
         $log->info("correct file: $correctFile out of the two: " . join(' | ', @files));
      }
      else {
         $log->error("$fileName could not find sequence aln file");
         die "did not find sequence alnfile\n";
      }
      my $alntxt = readAlnText($correctFile);
      push @{$result->[$i]}, $alntxt;
   }
   system("pwd");
   $log->info("QC result transformed");
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


=head2 readQC

 header:
   input fraction ends direction_ref direction_ass alnlen numgap_ref numgap_ass gaplen_ref gaplen_ass identity begin_ref end_ref begin_ass end_ass alntext); # add other columns

=cut
sub readQC {
   my $self = shift;

   my $file = $self->getQCFileName;

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
   close IN;
   $self->{qc} = \@results;
   print "QC read from $file\n";
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

####################################
############# result part ##########

=head2 uploadResult

 will upload result from result file into database. It also do some transformation
 to match the database schema. The original result file was designed for human
 consumption.

 return true if success false otherwise.

=cut
sub uploadResult {
   my $self = shift;
   
   $log->info("uploading mira result for " . $self->getProjectName . " to database...");
   my $workdir = $self->getWorkDirectory;
   chdir $workdir;
   eval {
      $self->readResult;
      $log->debug("result read from file");
      $self->transformResult;
      $log->debug("result transformed");
      if (!$self->uploadResultToDatabase) {
         return 0;
      }
   };
   if ($@) {
      $log->fatal("failed one of the steps: read result, transform result, or upload result to db");
      return 0;
   }
   else {
      chdir "-";
      $log->info("Reading result and loading into database done");
      return 1;
   }
}

=head2 uploadResultToDatabase 

 upload the result from file to database for the current project.
 return true of success false otherwise

=cut
sub uploadResultToDatabase {
   my $self = shift;

   my $result = $self->{result};
   my $dbh = $self->{dbh};
   my $request_id = $self->getRequestId;

   my $sth = $dbh->prepare("insert into result (input, fraction, input_coverage, contig_length, quality, coverage, head, tail, rawcut, request, sequence) values (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)");

   my $h = $result->[0];
   print "result header:\n", join("\t", @$h), "\n";
   eval {
      for (my $i=1; $i < @$result; ++$i) {
         $sth->execute(@{$result->[$i]});
      }
   };
   if ($@) {
      $log->fatal("Failed to upload result to database $@");
      $dbh->rollback or die $dbh->errstr;
      return 0;
   }
   else {
      $dbh->commit or die $dbh->errstr;
      $log->info("result loaded to database");
      return 1;
   }
}

sub transformResult {
   my $self = shift;

   my $result = $self->getResult;
   my $requestId = $self->getRequestId;
   my $resultDir = $self->getResultDirectory;
   my $seqreader = $self->getSequenceReader;
   $log->debug("transforming results ...");

   if (!$result) {
      die "got no result!\n";
   }
   $log->info(scalar(@$result) . " result rows");
   my $h = $result->[0];
   $log->debug("Old header:\n" . join(' | ', @$h));
   splice @$h, 8, 2;
   push @$h, ("request", "sequence");
   $log->debug("Result Header after header transformation:\n" . join("\t", @$h));
   $log->info(scalar(@$result) . " results to transform");

   for (my $i=1; $i < @$result; ++$i) {
      splice @{$result->[$i]}, 8, 2;
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
      push @{$result->[$i]}, ($requestId, $seq->sequence);
   }
   $log->info("result transformed for database upload");
}

sub getResult {
   my $self = shift;
   return $self->{result};
}

=head2 readResult

 read result into data structure [row data]
 The first row is the header.

=cut
sub readResult {
   my $self = shift;
   
   my $resultFile = $self->getResultfile;
   open IN, "<$resultFile" or die $!;
   my $line = <IN>; # header
   if (!$line) {
      $log->fatal("empty line\n");
      die "Empty line\n";
   }
   chomp $line;
   my @header = split /\t/, $line;
   if (@header != 11) {
      $log->fatal("headers:\n" . join(' | ', @header) . "\n" . scalar(@header)
         . " columns not 11 as expected!\n");
      $log->fatal("expecting the following columns: input frac  input_coverage length_largest quality  coverage head  tail  head_rc  tail_rc  rawcut");
      die;
   }

   my @results;
   push @results, \@header;
   while (<IN>) {
      chomp;
      my @row = split /\t/, $_;
      push @results, \@row;
   }
   $self->{result} = \@results;
   close IN;
   return \@results;
}


sub fraction2string {
   my $frac = shift;
   if ($frac == 1) {
      return 1;
   }
   return int(sprintf("%.4f", $frac)*10000);
}

1;
