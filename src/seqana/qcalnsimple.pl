#!/usr/bin/perl -w
#file: qcalnsimple.pl

# for picking out very poor alignment from simpe sequence matches

$infile = "cllow.aln";
$outfile = "highsimple.pair";

@badpat = qw(S{10,} T{10,} Q{10,} (T.A){9,} [TA]{15,} [SN]{20,} [SD]{20,} (T.{1,2}){10,} (G..){12,} (P.D){9,} (A+.{1,4}){20,} (SS.){10,} ((S+).{1,2}){10,} ([QE]+.){10,} ([TS]+.){10,} (Q+.{1,3}){9,} (TT.{4]){10,} (TTT.{3}){10,} (T.P){10,});
$i=0;
while ($ARGV[$i]) {
	if ($ARGV[$i]) { $infile = $ARGV[$i]; }
	$i++;
}

open IN, "<$infile" or die "Cannot open $infile $!\n";
open HIGH, ">$outfile" or die "cannot open $outfile for output\n";
open BADTAB, ">badpair.tab";  # bad
open GDTAB, ">goodpair.tab";  # good

readAln();
while ($_) {
	readAln();
}

####################################################


sub readAln { #read one alignment and return three arrays
	$_ = <IN>;
	while ($_ && !/^={9,}/) { $_ = <IN>; }
	if (eof) { return; }

	$pair = <IN>;
	$pair =~ s/\s+/\t/;
	while ($_ && !/^Score: /) { $_ = <IN>; }
	while ($_ ne "//\n") {
		$buff="";
		$_ = <IN>;
		$buff .= $_;
		/(\s+\d+ )([A-Z\-]+)(\s+\d+)/;
		$query = $2;
		$lnlen = length($query);
		$lead_space = length($1);
      #$tail_space = length($3);
		$_ = <IN>;  # read middle line
		$buff .= $_;
		$spacer = substr($_, $lead_space, $lnlen);
		$_ = <IN>;   # read target sequence
		$buff .= $_;
		$target = substr($_, $lead_space, $lnlen);
		$_ = <IN>;   # read empty line
		$buff .= $_;
		$_ = <IN>;
		while (!/length=/) {
			$buff .= $_;
			$query .= substr($_, $lead_space, $lnlen);
			$_ = <IN>;
			$buff .= $_;
			$spacer .= substr($_, $lead_space, $lnlen);
			$_ = <IN>;
			$buff .= $_;
			$target .= substr($_, $lead_space, $lnlen);
			$_ = <IN>; # read space
			$buff .= $_;
			$_ = <IN>;
		}
		$query =~ s/\s+\d+$//;
		$target =~ s/\s+\d+$//;
		$_ = <IN>;
		$_ = <IN>; #read align_len=223 nogap_align_len=202
		$identity_ln = <IN>;
		@arr = split / /, $identity_ln;
		$arr[0] =~ s/identity=//;
		$arr[1] =~ s/nogapIdentity=//;
		$arr[2] =~ s/simlarity=//;
		while ($_ ne "//\n") { $_ = <IN>; }
		$localIden=0;
		if (simple($query)  || simple($target)) 
		{
			#with high identity matches, not totally bad
			if ($arr[1] > 0.95 || highIdenSegSimple($spacer)) {
				print HIGH "$pair\n$query\n$spacer\n$target\n$identity_ln\n";
			}
			else { # bad
				#print "$pair\n$query\n$spacer\n$target\n$identity_ln\n";
				print "$pair\n$buff";
				print BADTAB $pair;
			}
		}
		elsif ($spacer =~ /\|{12,}/) {  # good
			#print GDTAB "$pair\n$query\n$spacer\n$target\n$identity_ln\n";
			print GDTAB "$pair\n$buff\n";
			#print GDTAB $pair;
		}
		elsif (highIdenSeg($spacer)) { #good
			print GDTAB "$pair\n$buff\n";
			#print GDTAB "$pair\n$query\n$spacer\n$target\n$identity_ln\n";
		}
		elsif ($localIden=highIdenWindow($spacer, 60, 0.6)) { #good
			print GDTAB "$pair local_identity=$localIden\n$buff\n";
		}
		else { # bad
			#print "$pair\n$query\n$spacer\n$target\n$identity_ln\n";
			print "$pair\n$buff";
			print BADTAB $pair;
		}
	}
}

sub simple {
	$str = shift;
	for ($i=0; $i<@badpat; $i++) {
		$pat = $badpat[$i];
		if ($str =~ /$pat/) { return 1; }
	}
	return 0;
}

sub highIdenSeg {
	my $sp = shift;
	if ( $sp =~ /\|{12,}/ 
		|| $sp =~ /\|{4,}.{1,20}?\|:\|{11,}.{1,2}?\|{1,}.{1,20}?\|{4,}/
		|| $sp =~ /\|{10,}.{1,200}?\|{10,}/ 
		|| $sp =~ /\|{7,}.{1,100}?\|{10,}.{1,100}?\|{7,}/ 
		|| $sp =~ /(\|{8,}.{1,200}?){3}/
		|| $sp =~ /(\|{6,}.{1,100}?){4}/
		|| $sp =~ /\|{9,}.{1,20}?\|{9,}:\|{9,}:\|{1,}?/
		|| $sp =~ /[\|:]{26,}.{1,10}?\|{4,}/
		|| $sp =~ /\|{6,}.{1,2}?\|{6,}:\|{1,}/
		|| $sp =~ /(\|{4,}.{1,100}?){6,}/
		|| $sp =~ /\|{4,}.{1,50}?(\|{6,}.{1,100}?){2}.{1,100}?\|{4,}/) {
		return 1;
	}
	return 0;
}
	
# for simple sequences should be higher strangency
sub highIdenSegSimple {
	my $sp = shift;
	if ( $sp =~ /\|{14,}/ 
			|| $sp =~ /\|{10,}.+\|{10,}/ 
			|| $sp =~ /(\|{8,}.+){3}/
			|| $sp =~ /(\|{3,}.){3,}.{1,10}?[\|:]{19,}/
			|| $sp =~ /(\|{6,}.+){4}/) {
		return 1;
	}
	return 0;
}

# (spacer_str, Window_size, Identity_cut)
# Identity_cut in decemal point format 0.01-0.99
sub highIdenWindow {
	my @arr = split //, $_[0];
	my $W = $_[1];
	my $idencut = $_[2];
	my $i=0;
	my $j=0;
	if (@arr < $W) { return 0; }
	my $sum = 0;
	for ($j = 0; $j<$W; $j++) { 
		if ($arr[$j] eq '|') { $sum++; } 
		elsif ($arr[$j] eq ':') { $sum  += 0.5; }
	} # compute initial value of sume of |
	$i++;
	while ($i+$W-1 < @arr) {
		if ($arr[$i-1] eq '|') { $sum--; }
		elsif ($arr[$i-1] eq ':') { $sum -= 0.5; }
		if ($arr[$i+$W-1] eq '|') { $sum++; }
		elsif ($arr[$i+$W-1] eq ':') { $sum += 0.5; }
		if ($sum/$W > $idencut) {
			return ($sum/$W);
		}
		$i++
	}
	return 0;
}
		
