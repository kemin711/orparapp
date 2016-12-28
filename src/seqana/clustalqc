### perl ###
##!/usr/bin/perl -w

# this program is for qc clustalw alignments produced by clustalw
# oter programs such as dialign also produce clustalw formats but they differ
#  and this program has to be modified


# I use a c program called qcaln to do the job, the perl script is
# just a front end program

@files = glob("*.aln");

$begin=0;
$end = @files;

$i=0;
while ($ARGV[$i]) {
	if ($ARGV[$i] eq "-b") { $begin = $ARGV[++$i]; }
	elsif ($ARGV[$i] eq "-e") { $end = $ARGV[++$i]; }
	$i++;
}

mkdir("good", 0777);
mkdir("bad", 0777);
mkdir("nocommoncol", 0777);
mkdir("goodpart", 0777);

open BAD, ">bad.log";
open GOOD, ">good.log";
open NOCOM, ">nocom.log";
open PART, ">goodpart.log;";
$goodCnt=0;
$badCnt=0;
$nocomCnt=0;
$goodpartCnt=0;

for ($i=$begin; $i<$end; $i++) {
#for ($i=0; $i<@files; $i++) {
	#print "$files[$i]: \n$quality\n";
	$stem = $files[$i];
	$stem =~ s/aln//;
	$diafile = $stem . "cw";
	$diaquality = `qcaln $diafile`;
	$cluquality = `qcaln $files[$i]`;   # clustal quality
	if ($cluquality =~ /VERY HIGH/ || $diaquality =~ /VERY HIGH/ || 
			$diaquality =~ /HIGH NOGAP/ || $cluquality =~ /HIGH NOGAP/) {
		system("cp $stem* good");
		print GOOD "$stem: \n$cluquality\n$diaquality\n";
		$goodCnt++;
	} 
	elsif ($cluquality =~ /HIGH PARTIAL/ || $diaquality =~ /HIGH PARTIAL/) {
		system("cp $stem* goodpart");
		print PART "$stem: \n$cluquality\n$diaquality\n";
		$goodpartCnt++;
	}
	elsif ($cluquality =~ /NOGAP COLUMN TOO FEW/ && $diaquality =~ /NOGAP COLUMN TOO FEW/) {
			system("cp $stem* nocommoncol");
			print NOCOM "$stem: \n$cluquality\n$diaquality\n";
			$nocomCnt++;
	}
	else {
			system("cp $stem* bad");
			print BAD "$files[$i]: \n$cluquality\n$diaquality\n";
			$badCnt++;
 	}
}
print "Good: $goodCnt  Partial Good: $goodpartCnt No common column: $nocomCnt  Bad: $badCnt\n";
