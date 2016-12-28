#!/usr/bin/perl -w
#this program converts genscan results into ace format

$infile = <>;
open RESULT "$infile" or die "$!";
#RESULT is the Genscan result

#open GS "predict.gs";  #contain the combined 
open ACE, ">genscan.ace";
	
	while (!/----/) { 
		$_ = <RESULT>;
	}
	$_ = <RESULT>; #read empty line
	$_ = <RESULT>; #read the first result line
	if (/^NO EXONS/) {
		print LOG "No predicted gene for $subseq\n";
	}
	else {  #only if there are some predictions
		while (!/^Predicted/) { #not reaching the end of gen.exon list
		################loop of different genes of the same seq########
			@LOL = ();  #store the exon boundaries
			$Init = 0;
			$Term = 0;
			$Sngl = 0;
			@plyA = ();
			@prom = ();
			$exonCount = 0;
			$exonLen = 0;

			while (/\w/) {  #not empty line
				@arr = split /\s+/;
				shift @arr; #removed the first empty field
				
				###PlyA and Prom are special cases
				if ($arr[1] eq "PlyA") { @plyA = ($arr[3], $arr[4]); }
				elsif ($arr[1] eq "Prom") { @prom = ($arr[3], $arr[4]); }
				else { #Intr or Init or Term exons
					$exonCount++;
					$exonLen += $arr[5];
					if ($exonCount == 1) { #first exon
						$strand = $arr[2];
						if ($strand eq "+") {
							$fr = $arr[6];
							$ph = $arr[7];
						}
						$gn = $arr[0];
						$gn =~ s/\..+//;
					}
					if ($strand eq "-") {
						$fr = $arr[6];
						$ph = $arr[7];
					}
					if ($arr[1] eq "Init") { $Init = 1; }
					elsif ($arr[1] eq "Term") { $Term = 1; }
					elsif ($arr[1] eq "Sngl") { $Sngl = 1; }
					push @LOL, [ ($arr[3], $arr[4]) ]; 
				}
				$_ = <RESULT>;
			}
			$_ = <RESULT>;
			while (/^\s+$/) {  #reading empty lines
				$_ = <RESULT>; 
			}

	################################output section#########################
			makeName(); #makes a key for the subsequence 
			$cs = 0; #codon_start

			if ($strand eq "+") {
				$begin = $LOL[0][0];
				$end = $LOL[$#LOL][1];
				$ss = $begin + 2;
				while (($ss+$cs)%3 != $fr) { $cs++; }
			}
			else {
				$begin = $LOL[$#LOL][0];
				$end = $LOL[0][1];
				$ss = $begin;
				while (($ss-$cs)%3 != $fr) { $cs++; }
			}
			$cs++;
			print ACE "\nSequence $mainseq\n";
			print ACE "Subsequence $subName  $begin  $end\n";
			if (@plyA) { print ACE "polyA_site @plyA \"GenScan predict\"\n"; }

			##############subsequence################
			print ACE "\nSequence $subName\n";
			print ACE "CDS_predicted_by GenScan\n";
			if ($strand eq "+") {
				for $ii (0..$#LOL) {
					$s = $LOL[$ii][0] - $begin + 1;
					$e = $LOL[$ii][1] - $begin + 1;
					print ACE "Source_Exons  $s  $e\n";
				}
			}
			else {
				for ($ii=$#LOL; $ii>=0; $ii--) {
					$s = $begin - $LOL[$ii][0] +1;
					$e = $begin - $LOL[$ii][1] + 1;
					print ACE "Source_Exons $s  $e\n";
				}
			}
			if ($Sngl || $Init) {
				print ACE "CDS 1 $exonLen\n";
			}
			else { 
				print ACE "Start_not_found\n"; 
				print ACE "CDS $cs  $exonLen\n";
				#cs is codon_start
			}
			if (!$Term && !$Sngl) { print ACE "End_not_found\n"; }
			
			print ACE "Method PredictCDS\n";
			#############Output section ends##############
		}
	
		$_ = <RESULT>;
		while (!/^>/) {  
		#get rid of empty lines between Predicted peptide sequence(s):
		#and the first peptide sequence
			$_ = <RESULT>;
		}
		while (!eof(RESULT) && /^>/) { #loop for different peptides
			$ii = index($_, "peptide_", 1) + 8;
			$tmp = substr($_, $ii);
			($gn, $aalen) = split /\|/, $tmp;
			$aalen =~ s/_aa//;
			makeName();
			$_ = <RESULT>;
			$pep = "";
			while (!eof(RESULT) && /[a-yA-Y]/) {
				$pep .= $_;
				$_ = <RESULT>;
			}
			$pep .= $_;  #the last line

			print ACE "\nSequence $subName\n";
			print ACE "Corresponding_protein  $subName\n";
			print ACE "\nProtein $subName\n";
			print ACE "Title \"GenScan prediction\"\n";
			if ($pep =~ /^X/) { print ACE "N_term_missing\n"; }
			if ($pep =~ /X$/) { print ACE "C_term_missing\n"; }
			print ACE "Peptide $subName $aalen\n";
			print ACE "\nPeptide $subName\n$pep\n";

			while (!eof(RESULT) && !/^>/) { 
				$_ = <RESULT>;
			}
		}
	}
	close RESULT;
	$i++;


sub makeName {
	$subName = $subseq . "." . $gn;
	$subName =~ s/all/gs/;
}


