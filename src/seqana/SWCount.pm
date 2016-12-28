package SWCount;

use Exporter;
use strict;
use Carp;

our @ISA=qw(Exporter);
our @EXPORT=qw(mergeGenotype mergeCodon mergeBase);

=head1 NAME

SWCount - Post Processing for the comparetoref program

=head1 DESCRIPTION

This package will host all postprocessing routines for the comparetoref program.

=head2  mergeGenotype

=cut
sub mergeGenotype {
   my $inpath=shift;
   my $oupath=shift;

   open IN, "<$inpath" or die $!;
   open OU, ">$oupath" or die $!;

   my %geno=();
   my $line=<IN>;
   print OU $line; # header
   $line=<IN>;
   while ($line) {
      chomp $line;
      my @row=split /\t/, $line;
      $geno{$row[0]} = $row[1];
      $line=<IN>;
   }
   foreach my $g (keys %geno) {
      if ($g =~ /[a-z]/) {
         if (exists $geno{uc($g)}) {
            $geno{uc($g)} += $geno{$g};
         }
         else {
            $geno{uc($g)} = $geno{$g};
         }
         delete $geno{$g};
      }
   }
   my $total=sumHash(\%geno);
   foreach my $g (sort { $geno{$b} <=> $geno{$a} } keys %geno) {
      print OU $g, "\t", $geno{$g}, "\t", sprintf("%.8f", $geno{$g}/$total), "\n";
   }
   close IN;
   close OU;
}

sub mergeCodon {
   my $inpath=shift;
   my $oupath=shift;

   open IN, "<$inpath" or die $!;
   open OU, ">$oupath" or die $!;
   my $line=<IN>;
   while ($line) {
      chomp $line;
      my @row=split /\t/, $line;
      my $newRow=mergec(\@row);
      print OU join("\t", @$newRow), "\n";
      $line=<IN>;
   }
   close IN;
   close OU;
}

sub mergeBase {
   my $inpath=shift;
   my $oupath=shift;

   open IN, "<$inpath" or die $!;
   open OU, ">$oupath" or die $!;
   my $line=<IN>;
   print OU $line;
   $line=<IN>;
   while ($line) {
      chomp $line;
      my @row=split /\t/, $line;
      my $newRow=mergeb(\@row);
      print OU join("\t", @$newRow), "\n";
      $line=<IN>;
   }
   close IN;
   close OU;
}

sub mergec {
   my $row=shift;
   my $pos=shift @$row;
   my %cc=();
   # parse one line
   for (my $i=0; $i<@$row; ++$i) {
      if ($row->[$i] =~ /^([A-Za-z-]{1,3}):(\d+) \(/) {
         $cc{$1}=$2;
      }
      else {
         die "abnormal codon: ", $row->[$i], "\n";
      }
   }
   # merge lower to upper
   foreach my $c (keys %cc) {
      if ($c =~ /[a-z]/) {
         if (exists $cc{uc($c)}) {
            $cc{uc($c)} += $cc{$c};
         }
         else {
            $cc{uc($c)} = $cc{$c};
         }
         delete $cc{$c};
      }
   }
   # get sum
   my $total=sumHash(\%cc);
   my @result=($pos);
   foreach my $c (sort { $cc{$b} <=> $cc{$a} } keys %cc) {
      my $frac=sprintf("%.8f", $cc{$c}/$total);
      push @result, "$c:" . $cc{$c} . "($frac)"; 
   }
   return \@result;
}

sub mergeb {
   my $row=shift;
   if (!hasLowerBase($row)) {
      return $row;
   }
   my $pos=shift @$row;
   my %bc=();
   for (my $i=0; $i<@$row; $i += 2) {
      if ($row->[$i+1] =~ /(\d+)\([0-9.e-]+\)/) {
         $bc{$row->[$i]}=$1;
      }
      else {
         die "improper base column: ", $row->[$i+1], "\n";
      }
   }
   foreach my $b (keys %bc) {
      if ($b =~ /^[a-z]$/) {
         if (exists $bc{uc($b)}) {
            $bc{uc($b)} += $bc{$b};
            delete $bc{$b};
         }
         else {
            print join(" | ", @$row), "\n";
            confess "$b has no corresponding uppler case for location $pos\n";
         }
      }
   }
   # get sum
   my $total=sumHash(\%bc);
   my @result=($pos);
   foreach my $n (sort { $bc{$b} <=> $bc{$a} } keys %bc) {
      my $frac=sprintf("%.8f", $bc{$n}/$total);
      push @result, $n, $bc{$n} . "($frac)"; 
   }
   return \@result;
}

sub sumHash {
   my $bc=shift;

   my $total=0;
   foreach my $b (keys %$bc) {
      $total += $bc->{$b};
   }
   return $total;
}

sub hasLowerBase {
   my $row=shift;
   for (my $i=1; $i<@$row; $i += 2) {
      if ($row->[$i] !~ /^[A-Z]$/) {
         return 1;
      }
   }
   return 0;
}

1;
