package UncoverPolyA;

use Exporter;
use strict;
use Peak;

our @ISA=qw(Exporter);
our @EXPORT=qw(generateComplexityArray findPeaks generateAPairArray maxPeakByArea cutPolyAPeakTail);

sub maxPeakByArea {
    my $peaks=shift; # reference to array of peaks
    my ($i, $maxi, $max);
    $max=-10000;
    for ($i=0; $i<@$peaks; $i++) {
        if ($peaks->[$i]->area > $max) {
            $maxi=$i;
            $max=$peaks->[$i]->area;
        }
    }
    if (wantarray) {
        return ($peaks->[$maxi], $maxi);
    }
    else {
        return $peaks->[$maxi];
    }
}

=head2 generateAPairArray

default word length 10

=cut
sub generateAPairArray {
    my $seqarr=shift;
    my $wd=shift;
    $wd=10 if (!$wd); # window sieze
    Peak->windowsize($wd);
    my @carray;
    my $L=@$seqarr-$wd-1;
    if ($L <= 0) {
        warn "sequence shorter than window length $wd\n";
        return $L;
    }
    my ($i,$j);
# initialize the change value
    my $change=0;
    for ($j=0; $j<$wd; $j++) { 
        if (isAAPair($seqarr->[$j], $seqarr->[$j+1])) {
            ++$change;
        }
        else { --$change; }
    }
    $i=0;
    while (1) {
        $carray[$i] =$change/2; # number from -5 to + 5
        if ($i>=$L) { last; }
        if (isAAPair($seqarr->[$i], $seqarr->[$i+1])) {
            --$change;
        }
        else { ++$change; }
        $j=$i+$wd;
        if (isAAPair($seqarr->[$j], $seqarr->[$j+1])) {
            ++$change;
        }
        else { --$change; }
        ++$i;
    }
    return \@carray;
}

###############################################
my (%config, $withpacnt, $DEBUG);
# return 1 for simple sequence
# return 2 for start bad
# return 3 for end bad
# return 5 for both start and end bad
# return 0 if sequence good
# return a negative number if too short
sub complexityTrim {
    my $strref=shift; # reference to sequence as string
    my $seqarray=shift; # sequence in array format
    my $wd=shift; # window sieze normally 20 or 40
    my $carray = computeComplexity($seqarray, $wd);

    if (length($$strref) <= $wd) {
        return (-1, $carray);
    }

    my $cstr=join('', @$carray); # complexity array and string
    if ($DEBUG) {
        print LOG "complexity string with window size: $wd\n";
        printSeq($cstr, *LOG);
        print LOG "\n";
    }
    #if ($cstr =~ /^[0-3]+$/) {
    if (isSimple($carray)) { # more accurate
        if ($DEBUG) { 
            print LOG ">>> > Garbage sequence < <<<<<<\n"; 
            printSeq($$strref, *LOG);
        }
        return (1, $carray);
    }

    my $rv=0;
    my $x=0;
#### processing the end of the sequence #####
# so that we can still use the start of the index
    my $junkstart=0;
    if ($cstr =~ /([01]+)$/g) {
        $junkstart=pos($cstr) - length($1);
    }
    elsif ($cstr =~ /([01]+)([0-6]+?)$/g) {
        if (length($1) / length($2) > $config{lhratio}
            || (length($1) > $config{llength} && length($2) < $config{hlength})) {
            $junkstart=pos($cstr) - length($1) - length($2);
            if ($junkstart < 10) {
                die "Junk end is too close to the start!\n";
            }
        }
    }
#### GGGGGGGOODSEQXJUNKSTART+++++
    if ($junkstart) {
        if ($DEBUG)  {
            if ($junkstart > length($$strref)) {
                print STDERR "$cstr\n";
                die "$junkstart ouside $$strref\n"
            }
            print LOG "Junk to be processed from the end $junkstart\n", 
                substr($$strref, $junkstart), "\n";
        }
        my @teststr=qw(AAAAAA CCCCCC GGGGGG TTTTTT N);
        my $minx=length($$strref);
        foreach my $ts (@teststr) {
            $x=index($$strref, $ts, $junkstart);
            if ($x > -1 && $x < $minx) { $minx=$x; }
        }
        if ($minx == length($$strref)) {
            if ($DEBUG) {
                print LOG "No Polymer or N in bad region. The sequence may not have a bad end after all.\n";
            }
        }
        else {
            if ($DEBUG) {
                print LOG "**** Bad end located at: ", $junkstart, "\n";
                print LOG "Junk sequence to remove:\n", substr($$strref, $minx), "\n";
            }
            $$strref=substr($$strref, 0, $minx);
            #print STDERR "cstr before\n$cstr\n";
            $cstr=substr($cstr, 0, $minx);
            #print STDERR "cstr after\n$cstr\n";
            $rv += 2;
        }
    }

####### working on begining of the sequence ##############
    my $junkend=0;
    if ($cstr =~ /^[01]+/g) {
        $junkend=pos($cstr); # index 1-after the junk region
    }
    elsif ($cstr =~ /^([0-5]+?)([01]+)/g) {
        if (length($2) / length($1) > $config{lhratio}
            || (length($2) > $config{llength} && length($1) < $config{hlength})) {
            $junkend=pos($cstr);
        }
    }
    if ($junkend) {
        if ($junkend > length($$strref)) {
            die "$junkend of string: \n", $$strref, "  \ncstr:\n$cstr\n";
        }
        my $tmp=substr($$strref, $junkend, $wd);
        if ($DEBUG) {
            print LOG "Junk region end index in complexity string $junkend\n";
            print LOG "last window of junk seq: $tmp\n";
        }

        while ($tmp =~ /N/g) { $x=pos $tmp; }
        while ($tmp =~ /A{6,}|C{6,}|G{6,}|T{6,}/g) {
            $x=pos($tmp);
        }
        if ($DEBUG) {
            if ($x) {
                print LOG "extended the junk in last window by $x\n";
            }
            else {
                print LOG "no extension of junk in the last window\n";
            }
        }
        my $endpos=$junkend + $x;
        if ($DEBUG) {
            print LOG "junk at start to throw away: ", substr($$strref, 0, $endpos), 
                "\nfrom\n$$strref\n";
        }
        if (length($$strref) <= $endpos) {
            die "endpos $endpos outside string\n$$strref\ncstr\n$cstr\n";
        }
        $$strref=substr($$strref, $endpos);
        $rv += 3;
    }
    return ($rv, $carray);
}

=head2 generateComplexityArray

    @param ($ntarray_ref, $wordsize)

    return a reference to the complexity array

=cut
sub generateComplexityArray {
    my $ntseq=shift; # reference to seq in array format
    my $wd=shift;
    #my @ntseq= split //, $$strref;
    my @carray;
    my $L=@$ntseq-$wd-1;
    if ($L <= 0) {
        #warn "sequence shorter than window length $wd\n";
        return undef;
    }
    my ($i,$j);
# initialize the change value
    my $change=0;
    for ($j=0; $j<$wd; $j++) { 
        if (differentBase($ntseq->[$j], $ntseq->[$j+1])) {
            ++$change;
        }
    }
    $i=0;
    while (1) {
# old method is n square, we will make it linear
        $carray[$i] =int(($change/$wd)*9);
        if ($i>=$L) { last; }
        if (differentBase($ntseq->[$i], $ntseq->[$i+1])) {
            --$change;
        }
        $j=$i+$wd;
        if (differentBase($ntseq->[$j], $ntseq->[$j+1])) {
            ++$change;
        }
        ++$i;
    }
    return \@carray;
}

=head2 findPeaks

    @param ($arrayRef, $base)

    $arrayRef is a reference to an array of numbers of integers
    $base is the start value of the peak

    -2-3-90123456998765432100-1-1-3
    the peak is from 123...21 if $base=0

=cut
sub findPeaks {
# base 5 max 9
    my $aref = shift; # array of numbers
    my $base=shift; # defult 0 or 5
    my $i=0;
    my @peaks=();
    my ($pk, $pkb, $area, $max, $maxi);
    while ($i<@$aref) {
        if ($aref->[$i] == $base) {
            while ($i < @$aref && $aref->[$i] == $base) { ++$i; }
            last if ($i >= @$aref);
            next if ($aref->[$i] < $base);
            $pkb=$i; # position after $base 
            $area=0; $max=0;
            while ($i < @$aref && $aref->[$i] >= $base) {
                $area += $aref->[$i] - $base;
                if ($aref->[$i] > $max) { 
                    $max=$aref->[$i]; $maxi=$i;
                }
                ++$i;
            }
            $pk = Peak->new($pkb, $i, $area, $max, $maxi);
            push @peaks, $pk;
        }
        else { ++$i; }
    }
    return \@peaks;
}

sub differentBase {
    my $b1=shift;
    my $b2=shift;
    if ($b1 ne $b2 && $b1 ne 'N' && $b2 ne 'N') {
        return 1;
    }
    else { return 0; }
}

sub simpleTrim {
    my $strref=shift;
    my $chg=0;
    while ($$strref =~ s/N[ACGT]{0,6}$//) {
        ++$chg; # changed
    }
    if ($$strref =~ /A{17,}[CGTN]$/) {
        chop $$strref;
        ++$chg; 
    }
    while ($$strref =~ s/^[ACGT]{0,6}N//) { ++$chg; }
    if ($$strref =~ /^[ACGN]T{17,}/) {
        $$strref = substr($$strref, 1);
        ++$chg;
    }
    if ($$strref =~ s/^[NT]{6,}// 
        || $$strref =~ s/^[ACGN]{1,2}[NT]{6,}//
        || $$strref =~ s/^[ACGT]{1,18}T{14,}//
        || $$strref =~ s/^[ACGT]{1,7}T{8,}//
        || $$strref =~ s/^[ACGT]{1,7}T{7,}CT{7,}//
     ) { ++$chg; }
    if ($$strref =~ s/[NA]{6,}$// || $$strref =~ s/[NA]{6,}[CGT]{1,2}$//) { ++$chg; }
    return $chg;
}

sub isSimple {
    my $carr_ref=shift;
    my %chash;
    my $low=0; 
    my $high=0; # low and high complexity region, by 4
    foreach my $c (@$carr_ref) {
        ++$chash{$c};
    }
    foreach my $c (keys %chash) {
        if ($c < 4) {
            $low += $chash{$c};
        }
        else {
            $high += $chash{$c};
        }
    }
    if ($low/($low+$high) > 0.9) {
        return 1;
    }
    return 0;
}

# helper for shortReadPolyA
sub isAAPair {
    my $b1=shift;
    my $b2=shift;
    if (($b1 eq 'A' || $b1 eq 'N') && ($b2 eq 'A' || $b2 eq 'N')) {
        return 1;
    }
    else { return 0; }
}

=head2 averageComplexity

    @param ($complexityArray, $b, $e)

    Where $b and $e a 0-based index of the range of the 
    form [$b, $e), $e is one pass the end of the range.

=cut
sub averageComplexity {
    my $comparr=shift;
    my $b=shift;
    my $e=shift;
    my $sum=0;
    if ($b > @$comparr) {
        #warn "range start $b befound end of complexity array\n";
        return 0;
    }
    if ($e > @$comparr || ! defined $e) {
        #warn "range end befound end of complexity array\n";
        $e=@$comparr;
    }
    return 0 if ($e==$b);
    if ($e < $b) {
        die "begin more than end\n";
    }
    for (my $i=$b; $i<$e; $i++) {
        $sum += $comparr->[$i];
    }
    $sum /= ($e-$b);
    return $sum;
}

=head2 cutPeak

    @param: $seqref, $seqarray_ref, $complexity_array_ref, PolyaPeak

    return true if the peak has been removed otherwise
    return false

=cut
sub cutPeak {
    my $seqref=shift;
    my $saref=shift;
    my $caref=shift;
    my $pk=shift;

    my $avgcx1 = averageComplexity($caref, 0, $pk->{begin}); 
# complexity array after peak length
    my $wd=PolyaPeak->windowsize();
    my $tlen = scalar(@$caref) - $pk->end - $wd;
    my $avgcx2 = averageComplexity($caref, $pk->end + $wd, scalar(@$caref));
    if ($pk->begin > $tlen
        && ($avgcx1 > $avgcx2-0.3 || 
            ($avgcx1 > 5 && ($tlen < 50 || $pk->begin > 150))
            )
    ) {
        my $cutsize=scalar(@$saref) - $pk->peakStart;
        $$seqref = substr($$seqref, 0, $pk->peakStart);
        $#{$saref} -= $cutsize;
        $#{$caref} -= $cutsize;
        return 1;
    }
    else {
        print LOG "subseq of peak + wd: ", 
            substr($$seqref, $pk->begin, $pk->width+PolyaPeak->windowsize),
                "\navg complexity before peak: $avgcx1",
                "\navg complexity after peak: $avgcx2\n";
        print LOG "Failed to cut away peak\n";
        $pk->show(*LOG);
        printSeq($$seqref, *LOG);
        print LOG join('', @$caref), "\n";

        ++$withpacnt;
        print SEQ ">seq$withpacnt\n";
        printSeq($$seqref, *SEQ);
        return 0;
    }
}
=head2 cutPolyAPeakTail

Remove the junk sequence after the polyA tail

-------------------| remove this part
CGCGAAAAAAAAAAAAAAACAAGCTGCAGTGT

=cut
sub cutPolyAPeakTail {
    my $seqref=shift;
    my $saref=shift;
    my $caref=shift;
    my $pk=shift;
    my $fh=shift;
    if (!$fh) { $fh=*STDOUT; }

    my $avgcx1 = averageComplexity($caref, 0, $pk->{begin}); 
# complexity array after peak length
    my $wd=Peak->windowsize();
    my $tlen = scalar(@$caref) - $pk->end - $wd;
    my $avgcx2 = averageComplexity($caref, $pk->end + $wd, scalar(@$caref));
    my $i;
    if ($pk->begin > 130 && ($avgcx1 > $avgcx2-0.1 || $avgcx1 > 5 )) {
        $i=$pk->peakStart;
        while ($i < @$saref && $saref->[$i] eq 'A') { $i++; }
        $$seqref = substr($$seqref, 0, $i);
        return 1;
    }
    else {
        print $fh "subseq of peak + wd: ", 
            substr($$seqref, $pk->begin, $pk->width+$wd),
                "\navg complexity before peak: $avgcx1",
                "\navg complexity after peak: $avgcx2\n";
        print $fh "Failed to cut away peak\n";
        $pk->show($fh);
        #printSeq($$seqref, *LOG);
        #print LOG join('', @$caref), "\n";

        #++$withpacnt;
        #print SEQ ">seq$withpacnt\n";
        #printSeq($$seqref, *SEQ);
        return 0;
    }
}

# input \$seq
sub shortReadPolyA {
    my $seqref=shift; # as reference
    my $seqarr=shift;
    my $comparr=shift; # complexity array
    my $wd=10; # window sieze
    if ($$seqref =~ s/(A{9,}[ACGTN]{1,9})$//) {
        my $polyalen = length($1);
        if ($$seqref !~ /A{9,}/) {
            return;
        }
        else {
            $#{$seqarr} -= $polyalen;
            $#{$comparr} -= $polyalen;
        }
    }

    my @carray;
    my $L=@$seqarr-$wd-1;
    if ($L <= 0) {
        warn "sequence shorter than window length $wd\n";
        return $L;
    }
    my ($i,$j);
# initialize the change value
    my $change=0;
    for ($j=0; $j<$wd; $j++) { 
        if (isAAPair($seqarr->[$j], $seqarr->[$j+1])) {
            ++$change;
        }
        else { --$change; }
    }
    $i=0;
    while (1) {
        $carray[$i] =$change/2; # number from -5 to + 5
        if ($i>=$L) { last; }
        if (isAAPair($seqarr->[$i], $seqarr->[$i+1])) {
            --$change;
        }
        else { ++$change; }
        $j=$i+$wd;
        if (isAAPair($seqarr->[$j], $seqarr->[$j+1])) {
            ++$change;
        }
        else { --$change; }
        ++$i;
    }
# find peak information
    $i=0;
    my @peaks=();
    my @discardedPeaks=();
    my $pk;
    while ($i<@carray) {
        if ($carray[$i] == 0) {
            my ($peakB, $peakE, $avgComplexity1, $avgComplexity2);
            while ($i < @carray && $carray[$i] == 0) { ++$i; }
            last if ($i >= @carray);
            next if ($carray[$i] < 0);
            $peakB=$i; # position after zero 
            my ($area, $max, $maxi);
            $area=0; $max=0;
            while ($i < @carray && $carray[$i] >= 0) {
                $area += $carray[$i];
                if ($carray[$i] > $max) { 
                    $max=$carray[$i]; 
                    $maxi=$i;
                }
                ++$i;
            }
            #$peakE=$i; # at zero
            $pk = PolyaPeak->new($peakB, $i, $area, $max, $maxi);
            if ($area > 25 || $pk->width > 7 || $pk->height > 1) {
                push @peaks, $pk;
            }
            elsif ($area > 4) {
                push @discardedPeaks, $pk;
            }
        }
        else { ++$i; }
    }
    if (@peaks == 0) { # no peaks saved but skipped
        if (@discardedPeaks > 0) {
            printSeq($$seqref, *LOG);
            print LOG join('', @carray), "\n";
            print LOG "no large peaks found\n",
                scalar(@discardedPeaks), " peaks discarded\n";
            foreach $pk (@discardedPeaks) {
                $pk->show(*LOG);
            }
        }
# no peaks detected no polyA
    }
    else {
        $i=$#peaks;
        while ($i > -1 && 
            cutPeak($seqref, $seqarr, $comparr, $peaks[$i]) )
        {
            --$i;
        }
        if ($i < $#peaks) { # actually removed peaks
            $#peaks = $i;
            if (@peaks > 0) {
# then use the lasrgest peak to curve
                my $maxAreaIdx=0;
                my $maxArea=0;
                for ($i=0; $i<@peaks; $i++) {
                    if ($peaks[$i] > $maxArea) {
                        $maxArea=$peaks[$i];
                        $maxAreaIdx=$i;
                    }
                }
                $pk=$peaks[$maxAreaIdx];
                print LOG "max peak after elimination:\n";
                $pk->show(*LOG); 
                print LOG "\nTrying to cut it again\n";
                cutPeak($seqref, $seqarr, $comparr, $peaks[$maxAreaIdx]);
            }
        }
    }
}

1;
