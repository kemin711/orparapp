package Fastaseq;

use strict;
use Carp;

my $width = 70; # line width used for printing
my %rctable = (
   A => 'T', C => 'G', G => 'C', T => 'A',
   R => 'Y', Y => 'R', S => 'S', W => 'W',
   K => 'M', M => 'K', B => 'V', D => 'H',
   H => 'D', V => 'B', N => 'N', '.' => '.',
   '-' => '-', '?' => '?',
   a => 't', c => 'g', g => 'c', t => 'a',
   r => 'y', y => 'r', s => 's', w => 'w',
   k => 'm', m => 'k', b => 'v', d => 'h',
   h => 'd', v => 'b', n => 'n'
   );

=head1 NAME 

Fastaseq - Fasta sequence simple operations.

=head1 DESCRIPTION

The basic operations for a fasta formated sequences and files.
It is simply represented by (name, description, sequence)
triplelet.  For complicated operation Bioseq should be used.

=head2 new

The constructor.

parameter($name, $description, $sequence)

 $sequence is the actual sequence string.

=cut
sub new {
   my $invocant = shift;
   my $name = shift;
   my $desc = shift;
   my $sequence = shift;
   my $class = ref $invocant || $invocant;
   my $self = { name => $name, description => $desc, sequence => $sequence };
   bless $self, $class;
   return $self;
}

sub clone {
   my $self = shift;
   if (!ref $self) {
      confess "clone method only works with objects of ", __PACKAGE__, "\n";
   }
   my $package = ref $self;
   return $package->new($self->name, $self->description, $self->sequence);
}


=head2 name

   getter and setter method.
   While argument is supplied it will be used to set the internal name of the sequence.

=cut
sub name {
   my $self = shift;
   if ($_[0]) {
      $self->{name} = $_[0];
   }
   return $self->{name};
}

=head2 description

getter and setter function for description.

=cut
sub description {
   my $self = shift;
   if ($_[0]) {
      $self->{description} = $_[0];
   }
   return $self->{description};
}

=head2 sequence

 getter and setter for sequence.
 sequence is simply a string.

=cut
sub sequence {
   my $self = shift;
   if ($_[0]) {
      $self->{sequence} = $_[0];
   }
   return $self->{sequence};
}

=head2 subseq

parameters($begin, $end)
If end is not given or negative then it is assumed the end of the sequence.

return a subsequence of the underlying sequence using 1-based index.


=cut
sub subseq {
   my $self = shift;
   my $b = shift;
   my $e = shift;
   my ($newName, $newDescription, $ss);
   if (defined $e && $e > 0) {
      $ss =  substr($self->{sequence}, $b-1, $e-$b+1);
      $newName = $self->name . '_' . $b . "_" . $e;
      $newDescription = $self->description . " subsequence from $b to $e";
   }
   else {
      $newName = $self->name . $b . "_" . "end";
      $newDescription = $self->description . " subsequence from $b to end\n";
      $ss = substr($self->{sequence}, $b-1);
   }
   return new Fastaseq($newName, $newDescription, $ss);

}

=head2 getFirst

parameter($length)

get the first $length of the parent seuqence 

return a new sequence.

=cut
sub getFirst {
   my $self = shift;
   my $len = shift;
   my $newName = $self->name . "first" . $len;
   my $newDescription = $self->description . " first " . $len;
   return new Fastaseq($newName, $newDescription, substr($self->{sequence}, 0, $len));
}

=head2 getLast

parameter($length)
return the last $length of sequence as a new sequence.

=cut
sub getLast {
   my $self = shift;
   my $len = shift;
   my $seq = $self->sequence;
   my $newName = $self->name . "last" . $len;
   my $newDescription = $self->description . " last $len";
   return new Fastaseq($newName, $newDescription, 
      substr($seq, length($seq) - $len, $len));
}

sub show {
   my $self = shift;
   print "name: $self->{name}\n",
      "description: $self->{description}\n",
      "sequence: $self->{sequence}\n";
}

=head2 revcomp 

Reverse complement this sequence.
The sequence will be changed, so do the name and description
to reflect the changes.

If the undelying sequence is not nucleic acid, the program will
die.

'U' in RNA should be converted to T for this function to work.
Otherwise, I will have to implement a derived class.

=cut
sub revcomp {
   my $self = shift;
   my @v = split //, $self->sequence;
   my $i=0; 
   my $j=$#v;
   while ($i < $j) {
      my $rc1 = $rctable{$v[$i]};
      my $rc2 = $rctable{$v[$j]};
      if (!defined $rc1 || !defined $rc2) {
         croak("Input sequence is not Nucleic Acids! " . $self->sequence . "\n");
      }
      $v[$i] = $rc2;
      $v[$j] = $rc1;
      ++$i; --$j;
   }
   if ($i == $j) {
      $v[$i] = $rctable{$v[$i]};
   }
   $self->{name} .= "rc";
   $self->{description} .= " reverse complement";
   $self->{sequence} = join ('', @v);
}

=head2 revcompCopy 

return a copy of the reverscomplement

=cut
sub revcompCopy {
   my $self = shift;
   my @v = split //, $self->sequence;
   my $i=0; 
   my $j=$#v;
   while ($i < $j) {
      my $rc1 = $rctable{$v[$i]};
      my $rc2 = $rctable{$v[$j]};
      $v[$i] = $rc2;
      $v[$j] = $rc1;
      if (!defined $rc1 || !defined $rc2) {
         die "Input sequence is not Nucleic Acids!\n";
      }
      ++$i; --$j;
   }
   if ($i == $j) {
      $v[$i] = $rctable{$v[$i]};
   }
   my $newName = $self->{name} . "rc";
   my $newDescription = $self->{description} . " reverse complement";
   return new Fastaseq($newName, $newDescription, join ('', @v));
}

=head2 relocateCutPoint

 return a new sequence with a new cut point.
 This operation only applies to circular DNA.

 parameter($position) is the new position 1-based index.
 The cut point will be starting at $position.
               pos 
 123....       |
 ==============------- old plasmid
 new one       -------==============

=cut
sub relocateCutPoint {
   my $self = shift;
   my $pos = shift;

   return new Fastaseq($self->name . "cut$pos",
              $self->description . " new cut point at $pos",
              substr($self->sequence, $pos - 1) . substr($self->sequence, 0, $pos - 1));
}

=head2 print

   write the sequence into a fasta foramted file handle.

=cut
sub print {
   my $self = shift;
   my $fh = shift;
   my $width = shift;
   if (!$width) { $width = 70; }

   print $fh ">", $self->name;
   if ($self->description) {
      print $fh " ", $self->description;
   }
   print $fh "\n";
   my $i = 0;
   while ($i < length($self->sequence)) {
      print $fh substr($self->sequence, $i, $width), "\n";
      $i += $width;
   }
}

=head2 printToFile

write the fasta formated sequence to a file.

=cut
sub printToFile {
   my $self = shift;
   my $fileName = shift;
   my $width = shift;
   if (!$width) {
      $width = 70;
   }
   open FOU, ">$fileName" or die "Cannot write to $fileName!\n";
   $self->print(\*FOU, $width);
   close FOU or die "Failed to close FOU. $!\n";
}

=head2 length

  return the length of the sequence.

=cut
sub length {
   my $self = shift;
   return length($self->sequence);
}

1;
