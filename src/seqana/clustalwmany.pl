### perl ###
##!/usr/bin/perl -w

# can run clustalw on many files in the same directory
# you have to start this program from the directory where the 
# *.pep files are located.
# due to a bug in clustalw, when single letter file x.pep is used 
# clustalw fails the file should be named C1,pep C2.pep ...
# this way the problem will be solved

# the result will be loaded into the database
# This is the database version, don't use this one
# Pg is the postgress perl interface
use Pg;
use dbinfo;

$infile_extn = "pep";  # default file extension file containing
                       # all sequences belonging to one cluster
$i= 0;
$b=0;  # begin index to align
$e=0;  # end index to align
$job="-align -endgaps ";  # default behavior of clustalw is to do alignment
##### clustalw cannot do tree and align at the same time
##### it must be run twice to get tree or align
$seqType="Protein"; # or DNA
$cluster_table = "pcluster";

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
		die "$ARGV[$i] not an option\n";
	}
	$i++;
}

# output file is named by clustalw automatically as infile.aln

$conn = get_orparadb;
$conn->exec("BEGIN");

@seqfiles = glob("*.$infile_extn");  #get all sequences in this dir
chomp @seqfiles;  # remove \n from each file name

if (!$e) { $e=@seqfiles; }

if ($MAKE_TREE) {
	$AllTreeFile = "ALL_in_one_file.bionj";
	open ALLTREE, ">$AllTreeFile";
}

for ($i = $b; $i < $e; $i++) {
	$seqfile = $seqfiles[$i];	
	#print "Doing Clustal $job on file $seqfile  # $i\n";
	$seqfile =~ /(\d+)\.$infile_extn/;
	$cluster_id = $1;

	$scorefile = nameOut($seqfile, "score");
	# $treeFile = nameOut($seqfiles[$i], "bionj");
	$dstFile = nameOut($seqfile, "dst");

	####### here is the major clustalw job ##############################
	$cmd = "clustalw -type=$seqType $job -INFILE=$seqfile > $scorefile";
	#print $cmd;
	system($cmd);
	######### clustalw has a bug when file x.pep is used as input and called with
	##### clustalw -aling -infile=x.pep it says that it could not read input file!

	#if ($job =~ /align/) { updateAlign(); }
	# update align model from the better of clustalw or dialign

	if ($MAKE_TREE) {           #need to modify bionj
 		# system("bionj $dstFile -o $treeFile");
		#### distance file *.dst is automatically generated if 
		# you call clustalw -infile=inf -tree -outputtree=dist
 		$tree = `bionj $dstFile 2>/dev/null`;
		$qresult = $conn->exec("UPDATE $cluster_table SET tree='$tree' WHERE id=$cluster_id");
		die $conn->errorMessage unless $qresult->resultStatus eq PGRES_COMMAND_OK;
		#print "\ntree result from bionj $tree\n";
		print ALLTREE "$cluster_id\n$tree\n";
 	}
}
$conn->exec("END");
print "\a\a\a Done \a\a\a!!!!!!!!!!!!!!!!!!!!!!!!!!\n";

sub nameOut 
{  # usage: nameOut(infile, ext_of_new_file);
   # returns new file with new extension, same stem
	my $ext=$_[1];
	my $newFile = $_[0];

	$newFile =~ s/$infile_extn$/$ext/;
	return $newFile;
}

# update the alignement in the database
# This should not be done because we will chosing the better align
# result from clustalw or dialign, only by tnen we need to add the 
# alignment to the database.
sub updateAlign {
	## produce multiple alignment db model
	print "cluster_id $cluster_id\n";
	$alnFile = nameOut($seqfiles[$i], "aln");
	($nameList, $consensus_len, $alignment) 
	                                 = split /\|/, `prodmodel $alnFile`;
	my $query;
	$query = "UPDATE $cluster_table SET consensus_len = $consensus_len, alignment = $alignment WHERE id = $cluster_id";
	$qresult = $conn->exec($query);
	die $conn->errorMessage unless $qresult->resultStatus eq PGRES_COMMAND_OK;
	#print "$nameList, $consensus_len \n $alignment\n";

	## update the order of the cds_id to be the same as those of gaplist
	$query = "UPDATE $cluster_table SET members = $nameList WHERE id = $cluster_id";
	$r = $conn->exec($query);
	die $conn->errorMessage unless $r->resultStatus eq PGRES_COMMAND_OK;
}

sub updateAlignment {
### simply insert the large alignment text into the database
### it may cause a lot of trouble
}
