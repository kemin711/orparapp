#!/usr/bin/perl -w

$nl=0;
$sumstart=0;
$sumstop=0;
$countstart=0;
$countstop=0;
while (<>) {
    #print $_;
    chomp;
    if ($nl > 5000) { last; }
    $idx=0;
    while (/TAG|TGA|TAA/g) {
        #print pos, " ";
        $d = pos($_) - $idx;
        $sumstop += $d;
        ++$countstop;
        if (!exists $dis{$d}) { $dis{$d} = 1; }
        else { ++$dis{$d}; }
        $idx=pos;
    }
    #print "\n";
    pos($_)=0;
    $idx=0;
    while (/ATG/g) {
        $d = pos($_) - $idx;
        $sumstart += $d;
        ++$countstart;
        if (!exists $start{$d}) { $start{$d} = 1; }
        else { ++$start{$d}; }
        $idx=pos;
    }
    ++$nl;
}

$meanstart=$sumstart/$countstart;
$meanstop=$sumstop/$countstop;

$sumstop=0;
open OU, ">stopdistance.tab";
my $k;
foreach $k (sort { $a <=> $b } keys %dis) {
    $sumstop += $dis{$k}*($k-$meanstop)*($k-$meanstop);
    print OU $k, "\t", $dis{$k}, "\n";
}
open OU, ">startdistance.tab";
$sumstart=0;
foreach $k (sort { $a <=> $b } keys %start) {
    $sumstart += $start{$k}*($k-$meanstart)*($k-$meanstart);
    print OU $k, "\t", $start{$k}, "\n";
}

print "mean start distance: $meanstart mean stop distance: $meanstop\n";
print sqrt($sumstart/$countstart), " as std start\n",
    sqrt($sumstop/$countstop), " as std stop\n";

