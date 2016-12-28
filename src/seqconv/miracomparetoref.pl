### perl ###

my $HOSTNAME;
BEGIN {
   $HOSTNAME = $ENV{HOSTNAME};
   if ($HOSTNAME eq 'foo') {
      push @INC, "/apps/sys/perl-5.18.1/lib/site_perl/5.18.1";
      push @INC, "/home/zhouke/foo/perlmod";
      push @INC, "/home/zhouke/foo/perlmod/lib/site_perl/5.18.1";
      $ENV{ORACLE_HOME} = "/u01/home/oracle/product/11.2.0/dbhome_1";
      $ENV{LD_LIBRARY_PATH}="/apps/sys/gcc-4.9.0/lib64:/apps/sys/lib";
   }
   else {
      push @INC, '/home/zhouke/perlmod';
      push @INC, '/usr/local/lib64/perl5/auto/DBD';
      push @INC, '/usr/lib/perl5/site_perl';
      push @INC, '.';
      $ENV{ORACLE_HOME} = "/data/kzoracle/product/12.1.0/dbhome_1";
      $ENV{LD_LIBRARY_PATH}="/data/kzoracle/product/12.1.0/dbhome_1/lib:/usr/local/lib:/usr/local/lib64:/usr/lib64:/usr/local/lib64/perl5/auto/DBD/Oracle";
      $ENV{PERL5LIB}="/usr/local/lib64/perl5/auto/DBD:/usr/lib/perl5/site_perl:/home/zhouke/perlmod:.";
   }
}

use strict;
use warnings;
use MiraOracleLoader;
use Log::Log4perl;

my $runmira;
if ($HOSTNAME eq 'bar') {
   $runmira = "/usr/local/orpara/bin/runmira";
}
else {
   $runmira = "/home/zhouke/foo/bin/runmira";
}
# the config file should be put into some central place in the production
# environment.
Log::Log4perl->init("/home/zhouke/src/proj/seqconv/log.conf");
my $log = Log::Log4perl->get_logger('miracomparetoref');
$log->info("starting to upload data to the miramill database");

my ($projname, $workDir);
my $i = 0;
while ($ARGV[$i]) {
   if ($ARGV[$i] eq "-w") { $workDir = $ARGV[++$i]; }
   elsif ($ARGV[$i] eq "--program") { $runmira = $ARGV[++$i]; }
   else {
      $projname = $ARGV[$i];
   }
   ++$i;
}
if (!$projname || !$workDir) {
   usage();
   die "you must give both project name and work directory name\n";
}

# 1 first we need to run the actual comparison part
$log->info("Doing comparison ...");
my $presentWorkingDirectory = `pwd`;
chomp $presentWorkingDirectory;
my $jobDirectory = "$workDir/$projname";
if ($presentWorkingDirectory ne $jobDirectory) {
   chdir $jobDirectory; 
}
$log->debug("Now in directory: $jobDirectory going to run $runmira\n");
my $inputFastq = $projname . ".fastq";
my $cmdstr = "$runmira -c $inputFastq";
$log->debug("Running $cmdstr");
system($cmdstr);
if ($?>>8) {
   $log->fatal("failed to run $cmdstr");
   die "Failed to run comparison part of mira QC\n";
}
else {
   $log->debug("$cmdstr done");
}
# 2 run the loading part
$log->info("compare done doing loading");
my $loader = MiraOracleLoader->new($projname);
# this is testing from europa on kraken file system
# so we have to be sure to use the right path
$loader->setRunDirectory($workDir);
if ($loader->upload) {
   $loader->markDone;
   $log->info("miracomparetoref done");
}
else {
   $log->fatal("miracomparetoref fail");
   die "Failed upload\n";
}

###################################################
sub getLogFileName {
   return "miracomparetoref.log";
}

sub usage {
   print STDERR "Usage miraload -w /net/kraken/ng10/kzassembly bar11\n",
      "  Options\n",
      "   --program path/to/runmira\n",
      "   argument project_name\n";
}

__END__

=head1 NAME

miracomparetoref - run the comparision part then to load assembly results to the database

=head1 DESCRIPTION

You first have to run the assembly pipeline, then run the comparison to reference part
with the runmira -c <projname>.fastq.  This will do the quality check.
Then you can run this program.


=head2 Options

 -w working directory. This is the default mira run place. Such as 
    /net/kraken/ng10/kzasembly, or /d1/work/assembly

 --program the absolute path for runmira.  Default is /usr/local/orpara/bin/runmira
   You can give a separate path. This is useful for production environment.

   The argument should be used for the project name, such as bar1, bar2, ...
