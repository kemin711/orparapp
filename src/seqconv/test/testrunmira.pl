### perl ### 
#This program will not be in the make system, but will be
#added to the hg system.

use File::Basename;
use Config::Simple;
use FastaReader;
use Fastaseq;

my $fileName = "bar1.fastq";

my ($name, $dir, $ext) = fileparse($fileName, '\.[a-z]+');

print "name: $name  dir: $dir  ext: $ext\n";

$fileName = "somefile";
($name, $dir, $ext) = fileparse($fileName, '\..*');
print "name: $name  dir: $dir  ext: $ext\n";


my $frac = 0.278973;
my $str = int(1000*sprintf("%.3f", $frac));
print $frac, " ", $str, "\n";
print sprintf("%.4f", 1)*10000, "\n";

collectAllLower('bar8_trim419cut.fastq', 0.765);
testConfig();
copyReferenceRC("bar3known.fas", "results");

sub testConfig {
   my $cfg = new Config::Simple("runmira.cfg");
   print $cfg->param('input'), "\n",
      " ascent", $cfg->param('ascent'), "\n",
      " ascend", $cfg->param('ascend'), "\n";
   if (!$cfg->param('ascend')) {
      print "ascend not in the config file\n";
   }
}

sub collectAllLower {
   my $infile = shift;
   my $fraction = shift; # top one 
   my $millage = int(1000*$fraction);

   my ($stem, $path, $suffix) = fileparse($infile, ".fastq");
   my $pattern = $stem . "_[1-9]*.fastq"; # for shell
   my @lowerFiles = `ls $pattern`;
   chomp @lowerFiles;
   #print join(" | ", @lowerFiles), "\n";
   $pattern = $stem . "_([0-9]{2,3}).fastq"; # for perl
   #print "pattern is $pattern\n";
   my @collection = ();
   foreach my $f (@lowerFiles) {
      if ($f =~ /$pattern/) {
         my $value = $1;
         if ($value < $millage) {
            #print "$f is good\n";
            push @collection, $f;
         }
         else {
            #print "$f is bad\n";
         }
      }
      else {
         die "unexpected file pattern $f\n";
      }
   }
   return \@collection;
}

sub copyReferenceRC {
   my $refpath = shift;
   my $dir = shift;

   # this should be used in the working directory
   my $reader = new FastaReader($refpath);
   my $reffas = $reader->next;
   my $outname = basename($refpath);
   $outname =~ s/\.fas/rc.fas/;
   $reffas->revcomp;
   $reffas->printToFile("$dir/$outname");
   print "rc reference file $outname written to directory $dir\n";
   return $outname;
}

