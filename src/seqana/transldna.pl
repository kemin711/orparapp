### perl ###
##!/usr/bin/perl -w

use Bioseq;

if (@ARGV<1) { 
	usage();
}
$REVERSE_COMPLEMENT = 0;
$frame = 1;
$i=0;
while ($ARGV[$i]) {
	if ($ARGV[$i] eq '-r') { $REVERSE_COMPLEMENT = 1 ; }
	elsif ($ARGV[$i] =~ /^\-(\d)$/) { $frame=$1; }
	elsif ($ARGV[$i] eq '-b') { $start = $ARGV[++$i]-1; }
	elsif ($ARGV[$i] eq '-e') { $finish = $ARGV[++$i]; }
	elsif ($ARGV[$i] eq '-f') { $FULL = 1; }
	elsif ($ARGV[$i] =~ /-[a-z]/) {
		die "This option ", $ARGV[$i], " Not legal in this program\n"
	}
	else { $infile = $ARGV[$i]; }
	$i++;
}
$dnaseq = getSeqFromFile($infile);
#print $dnaseq;
if ($finish) {
	#print "Range end specified\n";
	if ($start) {
		#print "Range start specified\n";
		$dnaseq = substr($dnaseq, $start, $finish-$start);
	}
	else { $dnaseq = substr($dnaseq, 0, $finish); }
}
else {
	if ($start) {
		$dnaseq = substr($dnaseq, $start);
	}
}

if ($REVERSE_COMPLEMENT) { 
	$rc_dnaseq = revcomp($dnaseq); 
	$pep = translate($rc_dnaseq, $frame);
}
else {
	$pep = translate($dnaseq, $frame);
}

if ($FULL) {
	printFull($dnaseq);
}
else {
	printSeq($pep, *STDOUT);
}

##################################################
sub usage {
	print "Usage: transldna infile [DNA inf fasta format]\n",
        "options: -b, -e, 1-based index inclusive. \n";
    exit(1);
}

sub printFull {
	my $tl1 = translate($dnaseq, 1);
	my $tl2 = translate($dnaseq, 2);
	my $tl3 = translate($dnaseq, 3);
	my $i=0;
	my @arr = split //, $dnaseq;
	my @arrt1 = split //, $tl1;
	my @arrt2 = split //, $tl2;
	my @arrt3 = split //, $tl3;
	my $it1 = 0;
	my $it2 = 0;
	my $it3 = 0;
	my $first_line=1;
	while ($i<@arr) {
# numbering
		for ($j=0; $j<66 && $i+$j<@arr; $j++) {
			if ($j % 20 == 0) { 
				$tmpstr = $i + $j + 1; # 1-based index
				$len = length($tmpstr);
				print $tmpstr; 
			}
			else { 
				if ($len < 2) {
					print ' '; 
				}
				else {
					--$len;
				}
			}
		}
		print "\n";

		for ($j=0; $j<66 && $i<@arr; $j++) {
			print $arr[$i++]
		}
		print " $i\n";
# frame 1
		for ($j=0; $j<22 && $it1<@arrt1; $j++) {
			print ' ', $arrt1[$it1++], ' ';
		}
		print "\n";
# frame 2
		#if ($line == 0) { print ' '; }
		print ' '; # start with second base
		for ($j=0; $j<22 && $it2<@arrt2; $j++) {
			#if ($j == 0) { print ' '; }
			print ' ', $arrt2[$it2++], ' ';
			#if ($j != 21) { print ' '; }
		}
		print "\n";
# frame 3
		if ($first_line) { 
			print '  '; 
			for ($j=0; $j<21 && $it3<@arrt3; $j++) {
				print ' ', $arrt3[$it3++], ' ';
			}
			$first_line = 0;
		}
		else {
			for ($j=0; $j<22 && $it3<@arrt3; $j++) {
				if ($j > 0) { print ' '; }
				print $arrt3[$it3++], ' ';
			}
		}
		print "\n\n";
	}
}

__END__

=head1 NAME

    transldna - translate DNA sequence into protein sequence.

=head1 DESCRIPTION

    This is a simple command line program to convert mRNA
    into peptide sequences.

=head1 OPTIONS and ARGUMENTS

    This program takes one argument as the input fasta file.

 -r reverse_complement the sequence before translation
 - <number> [1,2,3] to indicate the frame of translation
 -b <number> begin position of translation. 1-based index.
 -e <number> end of translation. Inclusive, 1-based index.
 -f flag for full translation
