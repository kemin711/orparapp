### perl ###
##!/usr/bin/perl -w 

# can run clustalw on many files in the same directory
# will not work on the database, take input from file
# should eliminate the ability to calculate bionj tree
# should just produce output of distance, The tree can be 
# easily produced separately

use Pg;

$CLUSTER_ID_IN_FILE_NAME = 0;

$infile_extn = "pep";  # default file extension file containing
                       # all sequences belonging to one cluster
$i= 0;
$b=0;   # start index of job, provide a sub range
$e=0;   # so that we can start several job at the same time
# inclusive 

# default behavior of clustalw is to do alignment
$job="-align -gapopen=2 -gapext=1 -pwgapopen=2 -pwgapext=1 -nopgap -nohgap ";  

$seqType="PROTEIN"; # PROTEIN or DNA  all cap, by clustalw
while ($ARGV[$i]) {
	if ($ARGV[$i] eq '-b') { $b = $ARGV[++$i]; }
  elsif ($ARGV[$i] eq '-x' || $ARGV[$i] eq '--fext') { 
    $infile_extn = $ARGV[++$i];  # input file extension
  }
	elsif ($ARGV[$i] eq '-e') { $e = $ARGV[++$i]+1; }
	elsif ($ARGV[$i] eq '-t') { 
		$job="-tree -tossgaps -kimura -outputtree=nj"; 
	}
	elsif ($ARGV[$i] eq '--type') { $seqType = $ARGV[++$i]; }
	elsif ($ARGV[$i] eq '-d') { 
    #to instruct clustaw to produce distance matrix
		$job="-tree -tossgaps -kimura -outputtree=dist"; 
		$MAKE_TREE = 1;  #bionj will make tree
	}
	else {
		print "$ARGV[$i] not an option\n";
	}
	$i++;
}

# align and distance calculation are two mutally exclusive 
# activities by the clustalw program

# output file is named by clustalw automatically as infile.aln
if ($CLUSTER_ID_IN_FILE_NAME) {
	use_clusterID();
}
else {
	print "Doing genearl clustalw with any input filename\n";
	use_anyseqfile();
}
###################################### 
# work in one directory
sub use_anyseqfile {
	@seqfiles = glob("*.$infile_extn");  #get all sequences in this dir
	if (!$e) { $e=@seqfiles; }

	$all_treefile = "ALL_TREE.bionj";
	open ALLTREE, ">$all_treefile";
	# This is for human to read
	# the following is for database input
	open BNJ, ">ALL_BIONJ.tab";

	for ($i = $b; $i < $e; $i++) {
		$seqfiles[$i] =~ /(.+)?\.$infile_extn$/;
		$fileStem = $1;
		$scorefile = $fileStem . ".score";
		$cmd = "clustalw -type=$seqType $job -INFILE=" . $seqfiles[$i] . " > $scorefile";
		#print STEDRR $cmd, "   ==> $i\n";
		if ($job =~ /^-tree/) {
			print STDERR "Making Tree for ";
		}
		else { print STDERR "Aligning "; }
		print STDERR $seqfiles[$i], "  ==> $i ...\n";
		system($cmd);

		if ($MAKE_TREE) {           #need to modify bionj
			#$dstFile = nameOut($seqfiles[$i], "dst");
			$dstFile = $fileStem . ".dst";
			# system("bionj $dstFile -o $treeFile");
			# clollect into individual files
			$tree = `bionj $dstFile 2>/dev/null`;
			if (!$tree) { 
				print "No tree produced for ", $seqfiles[$i], "\n";
				system("cat",  $seqfiles[$i]);
				die;
			}
			#$seqfiles[$i] =~ /(.+)?\.$infile_extn$/;
			print ALLTREE $fileStem, "\n$tree\n";
			print BNJ "$fileStem\t$tree";
			# the program bionj puts a line terminator 
		}
	}
	close ALLTREE;
	close BNJ;
}

sub use_clusterID {
	@cdsfiles = glob("*.$infile_extn");  #get all sequences in this dir
	if (!$e) { $e=@cdsfiles; }

	$AllTreeFile = "ALL_in_one_file.bionj";
	open ALLTREE, ">$AllTreeFile";
	for ($i = $b; $i < $e; $i++) {
		print "Doing Clustal $job on file $cdsfiles[$i]--$i\n";
		$cluster_id = substr($cdsfiles[$i], 1);
		$cluster_id =~ s/\..+//;
		$scorefile = nameOut($cdsfiles[$i], "score");
		# $treeFile = nameOut($cdsfiles[$i], "bionj");
		system("clustalw -type=$seqType $job -INFILE=$cdsfiles[$i] > $scorefile");

		if ($MAKE_TREE) {           #need to modify bionj
			$dstFile = nameOut($cdsfiles[$i], "dst");
			# system("bionj $dstFile -o $treeFile");
			$tree = `bionj $dstFile 2>/dev/null`;
			print ALLTREE "$cluster_id\n$tree\n";
		}
	}
}

sub nameOut 
{  # usage: nameOut(infile, ext_of_new_file);
   # returns new file with new extension, same stem
	my $ext=$_[1];
	my $newFile = $_[0];

	$newFile =~ s/$infile_extn$/$ext/;
	return $newFile;
}

