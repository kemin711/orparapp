### perl ###

# testing perl DBI
use strict;
use warnings;
use DBI;

my $datasource = 'dbi:Oracle:bioinfo';
my $dbh = DBI->connect($datasource, 'miramill', 'miramill', { RaiseError => 1 })
   or die "Connection to database $datasource failed: ", $DBI::errstr, "\n";
$dbh->{LongReadLen}=60372;
my $sth = $dbh->prepare("select id, name, description, refseq from request");
$sth->execute;
while (my @arr = $sth->fetchrow_array) {
   print join(' | ', @arr), "\n";
   #print $arr[3], "\n";
}
print "Test bioinfo done\n";

# test another database
# The following line does not work
#$datasource = 'dbi:Oracle:spandx@salus.pri.bms.com:1521';
#The following works
$datasource = 'dbi:Oracle:foo';
#                                     user     pass
my $dbh = DBI->connect($datasource, 'tcga', 'cancer', { RaiseError => 1 })
   or die "Connection to database $datasource failed: ", $DBI::errstr, "\n";
my $sth = $dbh->prepare("select id, acronym, definition from abbreviation");
$sth->execute;
while (my @arr = $sth->fetchrow_array) {
   print join(' | ', @arr), "\n";
   #print $arr[3], "\n";
}
print "test foo done\n";
