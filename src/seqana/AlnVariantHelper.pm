package AlnVariantHelper;

use strict;
use Bioseq;
use Carp;
use File::Basename;

our @ISA=qw(Exporter);
our @EXPORT=qw(extractVCFConsensus getFileStem);

=head1 NAME

AlnVariantHelper - helper routines for dealing with Sam and VCF files

=head1 DESCRIPTION

This is some of the itterface routines commonly used to deal with bam/sam
vcf/bcf files, but these routines are not implemented in the samtools and vcftools packages.

=head2 extractVCFConsensus

 input ($vcf_or_bcf_file, $output_fasta_file)

 the input file must ends with .vcf or .bcf
 the outputfile is not given will be created.

 return outputfilename this is the consensus.

=cut
sub extractVCFConsensus {
   my $file = shift;
   my $outfile = shift;
   if (!$outfile) {
      $outfile = substr($file, 0, index($file, '.')) . "_cons.fas";
      warn "no output file provided, we generated one for you $outfile\n";
   }

   ## open input source
   if ($file =~ /\.vcf$/) {
      open IN, "<$file" or die $!;
   }
   elsif ($file =~ /\.bcf$/) {
      open IN, "bcftools view $file |" or die $!;
   }
   else {
      confess "the input file must be *.vcf or *.bcf\n";
   }

   my $seq;
   while (<IN>) {
      next if (/^#/ || /^\s*$/); # skip comment or empty lines
      chomp;
      my @row=split /\t+/;
      $seq .= getConsensusBase(\@row);
   }
   close IN;

   open OU, ">$outfile" or die "Cannot open $outfile $!\n";

   my $prefix = substr($file, 0, index($file, '.'));
   my $consensusName = $prefix . '_consensus';
   print OU ">$consensusName\n";
   printSeq($seq, \*OU);
   print "consensus fasta sequence written to $outfile\n";
   close OU;
   return $outfile;
}

################# helper for extractVCFConsensus ########################
sub getConsensusBase {
   my $row = shift;

   my $AC=getAttribute($row->[7], "AC");
   my $REF=$row->[3];

   #print "$AC\n";
   my @acvalues = split /,/, $AC;
   my @altvalues = split /,/, $row->[4];
   if ($acvalues[0] == 2) {
      return $altvalues[0];
   }
   else {
      return $REF;
   }
}

sub getAttribute {
   my $str=shift;
   my $tag = shift;
   my @tmp = split /;/, $str;
   my %result = map { split /=/ } @tmp;
   return $result{$tag};
}

=head2 getFileStem

 return the base name, excluding suffix and path prefix

=cut
sub getFileStem {
   my $fpath=shift;
   my $sfx='';
   my $i = rindex($fpath, '.');
   #print "index of . in $fpath is $i\n";
   if ($i != -1) {
      $sfx = substr($fpath, $i);
   }
   #print "suffix is $sfx\n";
   my ($base, $path, $suffix)=fileparse($fpath, $sfx);
   return $base;
}

__END__

=head1 DESCRIPTION

This is a simple program to take vcf file and generate consensus.
If the input is compressed bcf then it should be able to detect it automatially
by first pipe the compressed form into the text format.



1;
