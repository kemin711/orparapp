### perl ###

use strict;
use warnings;
use MiraOracleLoader;
use Log::Log4perl;

# the config file should be put into some central place in the production
# environment.
Log::Log4perl->init("/home/zhouke/src/proj/seqconv/log.conf");
my $log = Log::Log4perl->get_logger('miraload');
$log->info("starting to upload data to the miramill database");


my ($projname, $workDir);
my $i = 0;
while ($ARGV[$i]) {
   if ($ARGV[$i] eq "-w") { $workDir = $ARGV[++$i]; }
   else {
      $projname = $ARGV[$i];
   }
   ++$i;
}
if (!$projname || !$workDir) {
   print "Usage miraload -w /net/kraken/ng10/kzassembly bar11\n";
   die "you must give both project name and work directory name\n";
}

my $loader = MiraOracleLoader->new($projname);
# this is testing from europa on kraken file system
# so we have to be sure to use the right path
$loader->setRunDirectory($workDir);
$loader->upload;

sub getLogFileName {
   return "miraloader.log";
}

__END__

=head1 NAME

miraload - to load assembly results to the database

=head1 DESCRIPTION

You first have to run the assembly pipeline, then run the comparison to reference part
with the runmira -c <projname>.fastq.  This will do the quality check.
Then you can run this program.


