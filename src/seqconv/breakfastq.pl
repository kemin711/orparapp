### perl ###

my $i=0;
my ($infile, $outfile, $numberOfPieces);
$numberOfPieces=2;
while ($ARGV[$i]) {
    if ($ARGV[$i] eq '-i') { $infile=$ARGV[++$i]; }
    elsif ($ARGV[$i] eq '-o') { $outfile=$ARGV[++$i]; }
    else {
        $numberOfPieces = $ARGV[$i];
    }
    ++$i;
}


breakup($infile, $numberOfPieces);

sub breakup {
    my ($infile, $numberOfPieces) = @_;

    my ($oustem, $suffix) = getFileHeadTail($infile);
    my $numseq = countNumberOfSequences($infile);

    open IN, "<$infile" or die $!;
    my $nsperfile = int($numseq/$numberOfPieces);
    print "one file has $nsperfile sequences\n";
    my @lines;
    my $i = 1;
    
    my $outfile = $oustem . '_' . int($i / $nsperfile) . '.' . $suffix;
    open OU, ">$outfile" or die $!;
    my $hasline = get4lines(\*IN, \@lines);
    while ($hasline) {
        print OU @lines;
        if ($i % $nsperfile == 0) {
            close OU;
            $outfile = $oustem . '_' . int($i / $nsperfile) . '.' . $suffix;
            open OU, ">$outfile" or die $!;
        }
        $hasline=get4lines(\*IN, \@lines);
        ++$i;
    }
}

sub get4lines {
    my ($fh, $arr) = @_;
    $arr->[0] = <$fh>;
    if (!$arr->[0]) {
        return 0;
    }
    $arr->[1] = <$fh>;
    $arr->[2] = <$fh>;
    $arr->[3] = <$fh>;
    return 1;
}

sub getFileHeadTail {
    my $file = shift;

    my $suffix='fastq';
    my $stem = $file;
    if ($file =~ /^(.+)\.(\w+?)$/) {
        $stem=$1;
        $suffix = $2;
    }
    return ($stem, $suffix);
}

sub countNumberOfSequences {
    my $file = shift;
    my $numline = `wc -l $file`;
    if ($numline =~ /^(\d+) .+/) {
        $numline = $1;
    }
    else { chomp $numline; }
    print $numline/4, " sequences\n";
    return $numline/4;
}


