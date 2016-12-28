### perl ###

use strict;
use warnings;
use MiraOracleLoader;
use Log::Log4perl;

Log::Log4perl->init("/home/zhouke/src/proj/seqconv/log.conf");
my $log = Log::Log4perl->get_logger('testmiraloader');
$log->info("starting to test the perl module for uploading data to the miramill database");

my $projname = $ARGV[0];
if (!$projname) {
   die "Usage testmiraloader.pl bar11\n";
}
my $loader = MiraOracleLoader->new($projname);

# this is testing from europa on kraken file system
# so we have to be sure to use the right path
$loader->setRunDirectory('/net/foo/ng10/kzassembly');

#$loader->fetchProjectInfo('bar10');
#$loader->showProjectInfo;
#$loader->uploadResult;
#if ($loader->uploadQC) {
#   $log->info("qc upload successful");
#}
#else {
#   $log->info("Failed to upload QC");
#}
$loader->upload;


sub getLogFileName {
   return "miraloader.log";
}
