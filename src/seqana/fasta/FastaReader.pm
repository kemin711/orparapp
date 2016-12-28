package FastaReader;

use strict;
use Fastaseq;
use Carp;

=head1 NAME

FastaReader - A simple fasta sequence file operator

=head1 DESCRIPTION

Do a few simple operations on fasta file.

=head2 new

Constructor($fastaFileName)

 usage: my $reader = new FastaReader($fas_file)
 my $fas_seq = $reader->next;

=cut
sub new {
   my $invocant = shift;
   my $fasfile = shift;

   my $class = ref $invocant || $invocant;
   my ($fh, $header, $fpos);
   # fhpos is the file handle position of reading the 
   # header of the first sequence.
   # This is used for looking for sequences in the entire file.
   if ($fasfile) {
      open ($fh, "<", $fasfile) or die $!, " cannot find file: $fasfile\n";
      $header = <$fh>;
      $fpos = tell $fh;
   }
   my $self = {fh => $fh, header => $header, fhpos => $fpos };
   bless $self, $class;
   return $self;
}

=head2 getFileHandle

return the file hadle connected to the fasta input file.

=cut
sub getFileHandle {
   my $self = shift;
   return $self->{fh};
}

=head2 useFile

Switch the reader to a new input file

=cut
sub useFile {
   my $self = shift;
   my $newFile = shift;
   my $fh;
   close $self->{fh} if $self->{fh};
   open $fh, "<$newFile" or croak($! . " Cannot open $newFile\n");
   $self->{fh} = $fh;
   $self->{header} = <$fh>;
}

=head2 next

return the next Fastaseq object.

=cut
sub next {
   my $self = shift;
   my $header = $self->{header};
   if (!$header) { return undef; }

   my $fh = $self->getFileHandle;
   my ($name, $desc, $seq);
   chomp $header;
   if ($header =~ /^>(.+)/) { # proper formated fasta file
      $header = $1; # getting rid of the first marker char '>'
      #print "header $header\n";
      if ($header =~ /(.+?)\s+(.+)/) {
         $name=$1;
         $desc=$2;
         #print "name: $name, desc: $desc\n";
      }
      else {
         $name = $header;
         $desc = "";
         #print "no description\n";
      }
      # read the actual sequences
      my $line = <$fh>;
      while ($line && $line !~ /^\s*$/ && $line !~ /^>/) {
         chomp $line;
         if ($line =~ /^\s*$/) {
            last;
         }
         if ($line !~ /^[a-zA-Z]+$/) {
            croak("fasta sequence residue not Bioseq character! $line\n");
         }
         $seq .= $line;
         $line = <$fh>;
      }
      if (!$seq) {
         confess("Empty input sequence!\n");
      }
      my $fastseq = new Fastaseq($name, $desc, $seq);
      $self->{header} = $line;
      return $fastseq;
   }
   else {
      confess("improperly formated fasta header: " . $header . "\n");
   }
}

=head2 getByName

get one Fastaseq object from the file.

 return undef if no fasta sequence in the file.

 Need to reset the file handle to the begeinning in case 
 the file pointer is in the middle of the file.


=cut
sub getByName {
   my $self = shift;
   my $seqname = shift;
   if (!$seqname) {
      confess("No name provided in getByName()!\n");
      die;
   }

   if (tell $self->{fh} > $self->{fhpos}) {
      seek $self->{fh}, $self->{fhpos}, 0;
   }

   while (my $fas = $self->next) {
      if ($fas->name eq $seqname) {
         return $fas;
      }
   }
   return undef;
}

=head2 getByMultipleNames

return multiple Fastaseq object in an array.
   If not a single one found, then the array is empty.
   Gives a warning if fewer sequences are found in the sequence file.

=cut
sub getByMultipleNames {
   my $self = shift;
   my $seqname = shift;
   my @result = ();
   # reset file handle to the start of first sequence
   # fhpos is the pos after reading the header of the 
   # first sequence.
   if (tell $self->{fh} > $self->{fhpos}) {
      seek $self->{fh}, $self->{fhpos}, 0;
   }
   while (my $fas = $self->next) {
      my $sn = $fas->name;
      if (grep { $_ eq $sn } @$seqname) {
         push @result, $fas;
      }
   }
   if (@result < @$seqname) {
      warn "only ", scalar(@result), " of ", scalar(@$seqname), " names found\n";
   }
   return \@result;
}

1;
