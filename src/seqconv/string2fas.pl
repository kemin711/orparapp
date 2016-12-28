### perl ###

# read a sequence file in pseudo fasta format and write a nice
# fasta file


my $infile=$ARGV[0];
open IN, "<$infile" or die $!;
open OU, ">nice.fas" or die $!;

my $ln = <IN>;
while ($ln) {
  if ($ln =~ /^>/) {
    chomp $ln;
    $ln =~ s/[\r\n]$//g;
    my $header=$ln;
    $ln = <IN>;
    my $seq;
    while ($ln && $ln !~ /^>/) {
      chomp $ln;
      $ln =~ s/[\r\n]$//g;
      $seq .= $ln;
      $ln = <IN>;
    }
    writeFasta($header, $seq, \*OU);
  }
  else {
    die "$ln should be starting with >\n";
  }
}

close IN;
close OU;


sub writeFasta {
  my $header=shift;
  my $seq = shift;
  my $fh = shift;
  print $fh $header, "\n";
  my $width=70;
  my $i=0;
  while ($i < length($seq)) {
    print $fh substr($seq, $i, $width), "\n";
    $i += $width;
  }
}

