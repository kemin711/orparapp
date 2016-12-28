package PolyaPeak;

=head2 new

    @param $begin, $end, $area, $hight, $firstHighestIdx
=cut
sub new {
    my $invocant=shift;
    my $class = ref($invocant) || $invocant;
    my $ref = { begin => $_[0], end => $_[1], area => $_[2], height => $_[3],
            firstHighest => $_[4]};
    bless $ref, $class;
    return $ref;
}

our $WINDOWSIZE=10;

sub windowsize { 
    if ($_[1]) { $WINDOWSIZE=$_[1]; }
    return $WINDOWSIZE; 
}

sub begin {
    my $self=shift;
    return $self->{begin};
}

sub end {
    my $self=shift;
    return $self->{end};
}

sub peakStart {
    my $self=shift;
    return $self->{firstHighest};
}

sub width {
    my $self=shift;
    return $self->{end}-$self->{begin};
}
sub height {
    my $self=shift;
    return $self->{height};
}

sub distance {
    my $self=shift;
    my $next=shift;
    return $next->{begin} - $self->{end};
}

sub show {
    my $self=shift;
    my $fh=shift;
    if (!$fh) { $fh=*STDERR; }
    print $fh "begin ", $self->{begin}, " end ", $self->{end},
        " area ", $self->{area}, " height ", $self->{height}, "\n";
}


1;
