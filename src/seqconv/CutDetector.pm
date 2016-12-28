package CutDetector;

#use Log::Log4perl;

my @transition = ( [0, 1, 2], 
                   [1, 1,-1],
                   [2, 1, 2] );


=head1 NAME

 CutDector - finite state automaton to detect the 
    cut point in the raw sequences

=head1 DESCRIPTION

 state: unknown or start 0, uncut 1, cut 2
 input: bad assembly 0, one large contig 1, two large contig 2

 State transition table:
           Input
 State | 0 | 1 | 2 | 
-------------------------
   0   | 0 | 1 | 2 |
------------------------
   1   | 1 | 1 | -1|
------------------------
   2   | 2 | 1 | 2 |
------------------------

=cut

sub new {
   my $invocant = shift;
   my $class = ref $invocant || $invocant;
   my $h = { state => 0 };
   bless $h, $class;
   return $h;
}

sub absorb {
   my $self = shift;
   my $input = shift;
   #my $log = Log::Log4perl->get_logger(ref $self);
   #$log->debug("state before absorb: " . $self->state);
   $self->state($transition[$self->state]->[$input]);
   #$log->debug("state after absorb: " . $self->state);
}

sub state {
   my $self = shift;
   if (defined $_[0]) {
      $self->{state} = $_[0];
   }
   else {
      return $self->{state};
   }
}

sub hasCut {
   my $self = shift;
   return $self->state == 2;
}

sub hasNoCut {
   my $self = shift;
   return $self->state == 1;
}

sub hasError {
   my $self = shift;
   return $self->state == -1;
}

1;
