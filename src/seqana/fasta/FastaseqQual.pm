package FastaseqQual;

use Fastaseq;
use strict;
#use Data::Dumper;
use Carp;

our @ISA = qw(Fastaseq);

=heqad1 NAME

FastaseqQual - fasta sequence with quality score.

=head2 new

 parameters($name, $description, $sequence, $quality)
 or
           ($fastaseq, $quality)

=cut
sub new {
   my $invocant = shift;
   my $self;
   if (@_ == 2) { # called with (Fastaseq, Quality)
      $self = $_[0]->SUPER::clone;
      bless $self;
      $self->quality($_[1]);
   }
   elsif ($_ == 4) {
      $self = $invocant->SUPER::new(@_); # made a copy of the current argument
      $self->{quality} = $_[3];
   }
   else {
      croak "expecting (name, description, sequence, ref_quality_array)\n or Fastaseq\n";
   }
   $self->{margin} = 12;
   $self->{flank} = 18;
   $self->{width} = 50;
   return $self;
}

sub quality {
   my $self = shift;
   if ($_[0]) {
      $self->{quality} = shift;
   }
   return $self->{quality};
}

sub margin {
   my $self = shift;
   if ($_[0]) {
      $self->{margin} = $_[0];
   }
   return $self->{margin};
}

sub flank {
   my $self = shift;
   if ($_[0]) {
      $self->{flank} = $_[0];
   }
   return $self->{flank};
}

sub width {
   my $self = shift;
   if ($_[0]) {
      $self->{width} = $_[0];
   }
   return $self->{width};
}

=head2 minQuality

 find the minimum quality from begining to end of the whole sequence.

=cut
sub minQuality {
   my $self = shift;

   if (!exists $self->{min}) {
      my $minq = 999999;
      my $minqi=0;
      my @Q;
      my $qq = $self->quality;
      for (my $i=0; $i < @$qq; ++$i) {
         if ($qq->[$i] < $minq) {
            $minq = $qq->[$i];
            $minqi = $i;
            @Q = ($i);
         }
         elsif ($qq->[$i] == $minq) {
            push @Q, $i;
         }
      }
      $self->{min} = $minq;
      $self->{mini} = \@Q;
   }
   if (wantarray) {
      return ($self->{min}, $self->{mini});
   }
   else {
      return $self->{min};
   }
}

sub minQualityInternal {
   my $self = shift;
   my $margin = shift;
   if ($margin) {
      $self->margin($margin);
   }
   else {
      $margin = $self->margin;
   }

   if (!exists $self->{internalMin}) {
      my $minq = 999999;
      my $minqi=0;
      my $qq = $self->quality;
      my @idx;
      for (my $i=$margin; $i < @$qq - $margin; ++$i) {
         if ($qq->[$i] < $minq) {
            $minq = $qq->[$i];
            $minqi = $i;
            @idx = ($i);
         }
         elsif ($qq->[$i] == $minq) {
            push @idx, $i;
         }
      }
      $self->{internalMin} = $minq;
      $self->{internalMini} = \@idx;
   }
   if (wantarray) {
      return ($self->{internalMin}, $self->{internalMini});
   }
   else {
      return $self->{internalMin};
   }
}

sub showInternalLowQuality {
   my $self = shift;
   my $margin = shift;
   my $fh = shift;
   if (!$fh) { $fh = \*STDOUT; }
   $self->margin($margin);
   print $fh "Internal low quality with margin: ", $self->margin, "\n";

   my ($minimumQuality, $minqIdx) = $self->minQualityInternal();
   $self->displayRegion($minimumQuality, $minqIdx, $fh);
}

sub showLowQuality {
   my $self = shift;
   my $fh = shift;
   if (!$fh) { $fh = \*STDOUT; }
   my ($minq, $minqIdx) = $self->minQuality;
   $self->displayRegion($minq, $minqIdx, $fh);
}

=cut displayRegion
 
 helper function to display the region for a given quality and
 all indices having such quality scores.

=cut
sub displayRegion {
   my $self = shift;
   my ($minimumQuality, $minqIdx) = @_;
   my $fh = $_[2];
   if (!$fh) { $fh = \*STDOUT; }

   print $fh "Flanking region(s) for quality $minimumQuality",
        ". Length of seq ", $self->length, "\n";
   my $flank = $self->flank;
   foreach my $ix (@$minqIdx) {
      print $fh "at 0-based index: $ix\n";
      if ($ix < $flank) {
         print $fh "index $ix inside the start flank region $flank, no left flank\n";
         print $fh " " x $ix, "| ", ($ix + 1), "\n";
         # subseq use 1-based index
         print $self->subseq(1, $ix + $flank + 1)->sequence, "\n";
      }
      elsif ($ix + $flank + 1 > $self->length) { # within the end 
         print $fh "index $ix inside the end flank region ",
            ($self->length - $flank), " no right flank\n";
         print $fh " " x $flank, "| ", ($flank + 1), "\n";
         # subseq use 1-based index
         print $fh $self->subseq($ix - $flank + 1)->sequence, "\n";
      }
      else { # ix fall inside the internal region excluding flank
         print $fh " " x $flank, "| ", ($flank + 1), "\n";
         # subseq use 1-based index
         print $fh $self->subseq($ix - $flank + 1, $ix + $flank + 1)->sequence, "\n";
      }
   }
}

sub showLowQualityBelow {
   my $self = shift;
   my $cutoff = shift;
   my $fh = shift;
   if (!$fh) { $fh = \*STDOUT; }

   print $fh "showing flanking regions for all qualities below $cutoff\n";
   my $lowq = $self->getQualityBelow($cutoff);
   foreach my $k (sort { $a <=> $b } keys %$lowq) {
      $self->displayRegion($k, $lowq->{$k}, $fh);
      print $fh "\n";
   }
}

sub getQualityBelow {
   my $self = shift;
   my $cutoff = shift;

   my $qq = $self->quality;
   my %lowq;
   for (my $i=0; $i < @$qq; ++$i) {
      if ($qq->[$i] < $cutoff) {
         push @{$lowq{$qq->[$i]}}, $i;
      }
   }
   return \%lowq;
}

=head2 printBoth

 A text-based pretty output of position, base, and score.
 Should get a better name.

=cut
sub show {
   my $self = shift;
   my $fh = shift;
   if (!$fh) { $fh = \*STDOUT; }

   my @s = split //, $self->sequence;
   my $q = $self->quality;
   my $i=0;
   my $w = $self->width;

   while ($i < @s) {
      my $j;
      my $ruler = ' ' x ($w*3);
      my $k = 0;
      for ($j=$i; ($j<@s && $j < $i + $w); ++$j) { # make ruler 
         if ($k % 10 == 0) {
            substr($ruler, $k*3, length($j)) = $j;
         }
         ++$k;
      }
      print $fh $ruler, "\n";
      for ($j=$i; ($j<@s && $j < $i + $w); ++$j) { # sequence character
         print $fh " ", $s[$j], " ";
      }
      print $fh "\n";
      for ($j=$i; ($j<@s && $j < $i + $w); ++$j) { # quality score
         print $fh sprintf("%2d", $q->[$j]), "|";
      }
      print $fh "\n";
      $i = $j;
   }
}

1;
