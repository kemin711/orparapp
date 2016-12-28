#!/usr/local/bin/perl -w 

#mainly to collect the annotation results from
#gtf to ace.  GTF files are scattered all over the place
#it's a whole mess

#$annot = "/data/data5/cg/h/annot.5";
$annot = "/data/data5/cg/h/affy.3";

#@chrDir = (1 .. 22, UL, X, Y);
@chrDir = (1 .. 22, X, Y);  #the UL directory has no *.mrg.gtf
#@chrDir = (1); #for testing  #the UL directory has no *.mrg.gtf

open ACE, ">annot.ace";
open ACESUB, ">annotSub.ace";

foreach $chr (@chrDir) {
	$dir = $annot . "/" . $chr;
	#@ntdir = <$dir/gb:NT_*>;  # for NT_ complete annot.5 files
	@seqdir = <$dir/ctg*>;

	foreach $subdir (@seqdir) {
		@arr = split /\//, $subdir;
		$fileName = $arr[7];
		my $sn = $fileName;
		my $DNAfile =$fileName . ".fa";
		my $PEPfile = $fileName . ".mrg.aa.fa";
		my $originalFile = $subdir . "\/" . $DNAfile;
		system("cp $originalFile $DNAfile");
		$originalFile = $subdir . "\/" . $PEPfile;
		system("cp $originalFile $PEPfile");
		$fileName .= ".mrg.gtf";
		$oldFile = $subdir . "\/" . $fileName;
		open IN, "<$oldFile" or die "cannot open $oldFile\n";
		#print "\nReading file $sn ------------>\n";
      convert2ace(); 
		close(IN);
	}
}
print "Done successfully\n";

sub convert2ace { #each file each sequence assumed
	my $ln = <IN>;
	while ($ln && $ln !~ /gb:NT_/) { $ln = <IN>; }
	if (!$ln) { 
		#print "No annotation!!\n";
		return; 
	}
	$ln =~ /^gb:(NT_\d+)/;
	$seqName = $1;
	#print "\nProcessing Sequence $seqName...\n";
	print ACE "\n\nSequence \"$seqName\"\n";

	my $CDSCnt = 0; 
	my $exonCnt = 0;

	#%cdsObj = ();
	#$cdsObj{"gene_id"} = "";  #starting condition
	#because the order of start/stop codon in gtf is
	#not assured to occur before the CDS listing is over
	#the output of the CDS should be deffered to the end of the 
	#input

	@cds=();  #to hold all CDS objects
	$rcds = {};
	$rcds->{"gene_id"} = "";  #starting condition
	$rcds->{"type"} = "CDS";
	%mRNA = ();
	$mRNA{"gene_id"} = "";
	$mRNA{"type"} = "mRNA";
	@startList = (), @stopList = ();
	%gene = ();  #to collect all gene_id gene_name pair

	while ($ln && $ln =~ /^gb:NT_/) {
		$ln =~ s/gb://g;
		my @arr = split /\t/, $ln;
		shift @arr;
		shift @arr;
		chomp($arr[6]);
		my %attrHash = ();
		my @attrArr = split /\"; /, $arr[6];
		my	$xx = shift @attrArr;
		while ($xx) {
			if ($xx =~ /exon_number (\d+)/) {
				$attrHash{"exon_number"} = $1;
			}
			else {
				$xx =~ /^(.+?) \"(.+)/;
				if (!$2) {
					print "format error:\n$ln\n$xx\n";
					goto NEXTLINE;
				}
				#$yy = $2;
				#if ($xx =~ /protein_id/) {
				#	$yy =~ s/.\d+//;
				#}
				#$attrHash{$1} = $yy;
				$attrHash{$1} = $2;
			}
			$xx = shift @attrArr;
		}
		if ($arr[0] eq "CDS") {
			$CDSCnt++;
			if (!($gid = $attrHash{'gene_id'})) {
				$attrHash{'gene_id'} = $seqName . $attrHash{'gene_name'};
			}  #making up an gene_id
				
			if ( $rcds->{"gene_id"} ne $attrHash{"gene_id"}) {
				if ($rcds->{"gene_id"} ne "") {  #new object 
					if ($rcds->{"strand"} eq "+") {
						$rcds->{"total_exons"} = $rcds->{"exon_number"} + 1;
					}
					if ($rcds->{'strand'} eq "-") {
						$rcds->{'frame'} = $arr[5];
					}
					push(@cds, $rcds);    $rcds = {};
				}
				loadHash(\%attrHash, \@arr, $rcds);
				if ($rcds->{'strand'} eq "+") {
					$rcds->{"frame"} = $arr[5];
				}  #the first coding exon contains the frame
				#in the plus strand
			}
			else {
				push( @{$rcds->{"source_exon"}}, [$arr[1], $arr[2]] );
				$rcds->{"exon_number"} = $attrHash{"exon_number"};
			}
		}  #CDS if
		elsif ($arr[0] eq "exon") {
			$exonCnt++;
			if ( $mRNA{"gene_id"} ne $attrHash{"gene_id"}) {
				if ($mRNA{"gene_id"} ne "") { 
					if ($mRNA{"strand"} eq "+") {
						$mRNA{"total_exons"} = $mRNA{"exon_number"} + 1;
					}
					writemRNAToAce();  %mRNA = ();
				}
				loadHash(\%attrHash, \@arr, \%mRNA);
			}
			else {
				push( @{$mRNA{"source_exon"}}, [$arr[1], $arr[2]] );
				$mRNA{"exon_number"} = $attrHash{"exon_number"};
			}
		}  #mRNA if
		elsif ($arr[0] eq "start_codon") {
			#print "loading start codon $ln\n";
			loadCodon(\@startList, \%attrHash, \@arr);
		}
		elsif ($arr[0] eq "stop_codon") {
			loadCodon(\@stopList, \%attrHash, \@arr);
		}
		elsif ($arr[0] eq "intron" || $arr[0] eq "gap"
		|| $arr[0] eq "CDS_insert") {
			#intron not considered, redundant information with exon
			#print "$arr[0] $arr[6]\nNot being considered\a\n";
		}
		else {
			print "$arr[0] $arr[1] $arr[6]\n new feature\n";
			die;
		}
NEXTLINE:		$ln = <IN>;
	}
	if ($CDSCnt>0) {
		push(@cds, $rcds);
		writeAllCDS();
	}
	if ($exonCnt>0) {
		writemRNAToAce();
	}
	foreach $gid (keys(%gene)) {
		print ACE "\nPGene $gid\n";
		my $gn = $gene{$gid};
		print ACE "Name \"$gene{$gid}\"\n" if ($gn ne "");
	}
}

#loadHash(\%attrHash, \@arr, \%mRNA);
#or                          \%CDS
sub loadHash {
	my $rh = $_[2];
	my $rfeatArr = $_[1];
	%$rh = %{$_[0]}; #direct assignment
	$rh->{"strand"} = $rfeatArr->[4];
	$rh->{"source_exon"}->[0] = [$rfeatArr->[1], $rfeatArr->[2]];
	if ($rh->{"strand"} eq "-") {
		$rh->{"total_exons"} = $rh->{"exon_number"} + 1;
	}
	if ($rfeatArr->[0] eq "CDS") {
		$rh->{'type'} = "CDS";
	}
	elsif ($rfeatArr->[0] eq "exon") {
		$rh->{'type'} = "mRNA";
	}
#	print "protein id $rh->{'protein_id'} before\n";
#	$rh->{"protein_id"} =~ s/\.\d+//;
#	print "protein id $rh->{'protein_id'} after\n";
}

sub writeAllCDS {
	foreach $c (@cds) {
		#print "Processing CDS gene_id= $c->{'gene_id'}\n";
		if (!$c->{'gene_id'}) {
			die "Sequence $seqName got no gene_id\n";
		}

		my $exonCnt = @{$c->{"source_exon"}};
		my $objStart = $c->{"source_exon"}->[0]->[0];
		my $objEnd = $c->{"source_exon"}->[$exonCnt-1]->[1];
		my $subSeqName = $c->{"protein_id"};
		if ($subSeqName) {
			$subSeqName =~ s/\.\d+//; #removing .2 
		}
		else {
			print  "Sequence $seqName have no protein_id, PSEUDO?\n";
			return;
		}
		#$c->{"protein_id"} = $subSeqName;
		#using protein_id as sequence name
		my $cdslen = 0;
		my ($i, $f, $s);

		print ACESUB "\nSequence $subSeqName\nMethod cdsobj\n";
		printFeat(*ACESUB, $c); 
		if ($c->{"strand"} eq "+") {
			print ACE "Subsequence $subSeqName $objStart $objEnd\n";
			foreach $ss ( @{$c->{"source_exon"}} ) {
				($xx, $yy) = @{$ss};
				$s = $xx - $objStart + 1;
				$f = $yy - $objStart + 1;
				$cdslen += ($f - $s + 1);
				print ACESUB "source_exons $s $f\n";
			}
		}
		else {
			print ACE "Subsequence $subSeqName $objEnd $objStart\n";
			for ($i = $exonCnt - 1; $i >=0; $i--) {
				($xx, $yy) =  @{$c->{"source_exon"}->[$i]};
				$s = $objEnd - $yy + 1;
				$f = $objEnd - $xx + 1;
				$cdslen += $f - $s + 1;
				print ACESUB "source_exons $s $f\n";
			}
		}
		my $startFound = findCodon(\@startList, $c);
		my $stopFound = findCodon(\@stopList, $c);

		if ($startFound eq "FALSE") {
			#die "no start_codon found\n";
			print ACESUB "Start_not_found\n";
			#my $fr = $c->{"frame"};
			my $CDS_start = $c->{"frame"} + 1;
			#my $CDS_start = $fr + 1;
			print ACESUB "CDS $CDS_start $cdslen\n";  
		}
		else {
			print ACESUB "CDS 1 $cdslen\n";  #not partial 
		}
		if ($stopFound eq "FALSE") {
			print ACESUB "End_not_found\n";
		}
	}
}


sub writemRNAToAce {
	#print "Processing mRNA gene_id= $mRNA{'gene_id'}\n";
	my $exonCnt = @{$mRNA{"source_exon"}};
	my $objStart = $mRNA{"source_exon"}->[0]->[0];
	my $objEnd = $mRNA{"source_exon"}->[$exonCnt-1]->[1];
	my $subSeqName = $mRNA{"transcript_id"};
	if (!$subSeqName) {
		print "mRNA has no transcript_id in $seqName\n";
		$subSeqName = $mRNA{'gene_id'};
	}
	$subSeqName =~ /(.*?NT_\d+?)\.\d+(-.+)/;
	$subSeqName = $1 . $2;
	my $mRNAlen = 0;
	my ($i, $s, $f);

	print ACESUB "\nSequence $subSeqName\nmRNA\nMethod mRNAobj\n";
   printFeat(*ACESUB, \%mRNA); 

	if ($mRNA{"strand"} eq "+") {
		print ACE "Subsequence $subSeqName $objStart $objEnd\n";
		foreach $ss ( @{$mRNA{"source_exon"}} ) {
			($xx, $yy) = @{$ss};
			$s = $xx - $objStart + 1;
			$f = $yy - $objStart + 1;
			$mRNAlen += ($f - $s + 1);
			print ACESUB "source_exons $s $f\n";
		}
	}
	else {
		print ACE "Subsequence $subSeqName $objEnd $objStart\n";
		for ($i = $exonCnt - 1; $i >=0; $i--) {
			($xx, $yy) =  @{$mRNA{"source_exon"}->[$i]};
			$s = $objEnd - $yy + 1;
			$f = $objEnd - $xx + 1;
			$mRNAlen += $f - $s + 1;
			print ACESUB "source_exons $s $f\n";
		}
	}
}

sub loadCodon {
#loadCodon(\@stopList, \%attrHash, \@arr);
	my $ra = $_[0]; #reference to codon array (start or stop)
	my $rarr = $_[2];
	my $rh = {};
	%$rh = %{$_[1]};
	$rh->{"start"} = $rarr->[1];
	$rh->{"end"} = $rarr->[2];
	$rh->{"strand"} = $rarr->[4];
	push(@$ra, $rh);
}

sub printLoH {
	my $ra = $_[0];
	foreach $st (@$ra) {
		print "*****************\n";
		foreach $key (keys(%$st)) {
			print "$key --> $st->{$key}\n";
		}
	}
}

#$startFound = findCodon(\@startList, $ref_to_CDS-Obj);
sub findCodon {
	my $ra = $_[0];
	my $rc = $_[1];
	my $i;
	if (@$ra < 1) {
		print "no entry in startList or stopList\n";
		return "FALSE";
	}

	for ($i=0; $i < @$ra; $i++) {
		if ($rc->{"protein_id"} eq $ra->[$i]->{"protein_id"}) {
			return "TRUE";
		}
	}
	return "FALSE";
}

sub printFeat {
#printFeat(*ACESUB, \mRNA or \CDS object) {
	local (*F) = $_[0];
	my $rh = $_[1]; #the hash object
	my ($gn, $gi);

	$rh->{'gene_id'} =~ /(.*?NT_\d+?)\.\d+(-.+)/;
	$rh->{'gene_id'} =  $1 . $2;

	if ($gi = $rh->{'gene_id'}) {
		if ($gn = $rh->{'gene_name'}) {
			$gene{$gi} = $gn;
		}
		else {
			print "Sequence $seqName gene_id $gi no gene_name******\n";
			$gene{$gi} = "";
		}
	}
	else {
		print "No gene_id***********\n"; 
	}
	my $p;
	if ($p=$rh->{'product'}) {
		print F "Title  \"$rh->{'product'}\"\n";
	}
	else {
		print "Sequence $seqName:\nprotein_id $rh->{'protein_id'} got no product****\n";
	}
	my $prt = $rh->{'protein_id'};
	if ($prt) {
		$prt =~ s/\.\d+//;
		#print F "Corresponding_protein  \"$rh->{'protein_id'}\"\n";
		print F "Corresponding_protein  \"$prt\"\n";
	}
	print F "$rh->{'type'}of_Gene \"$rh->{'gene_id'}\"\n";
}

