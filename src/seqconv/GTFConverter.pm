package GTFConverter;

use strict;

=head1 NAME

GTFConverter - convert GTF into Various formats mainly table

=head1 DESCRIPTION

The GTF format from Broad Institue is extremely simple
contain only exon, CDS -> transcript -> gene mapping
with start and stop codon information.

No extra information was added.

I will use this perl module to convert GTF formats into mainly
relational database schema which in turn could be easily changed
to other formats.

This is a copy and paste of GFF3Converter.  In the futrue,
this should be combined with GFF3Converter.

I will output table text formats no database interaction.

Using NCBI as test case 1

There could be more than one organism in one distribution
For example:

genomicid   strand  start   end mol_type    organism    chromosome  attribute
NC_008368.1 +   1   56814   genomic DNA Ustilago maydis \N  db_xref=taxon:5270; organelle=mitochondrion
NW_100966.1 +   1   30180   genomic DNA Ustilago maydis 521 1   db_xref=taxon:237631; strain=521

The download has not information as to how to construct the
sequence into contigs and contigs into scaffole and scaffold
into chromosomes!

the first test case is from broad institue.
Their ugly undocumented format.

The Broad GTF is a lot simpler than the GFF3 of NCBI

This module is written particularly for one database from Broad.

=cut

my %HTMLENCODING=('%20' => " ", '%2C' => ',', '%5B' => '[', '%5D' => ']'
    );

=head2 new

    usage my $gf3formater = GFF3Converter->new;

    This method will create several filehandles based on the
    attribute column 

    These handles will stay open during the object's life.
    So you can process multiple files and save then into
    the same streams.

=cut
sub new {
    my $invocant=shift;
    my $class=ref $invocant || $invocant;
    my $ref={
        column => { 
            CDS => ['genomicid', 'strand', 'phase', 'start', 'end', 'transcript_id'],
            exon => ['genomicid', 'strand', 'start', 'end', 'transcript_id'],
            start_codon => ['genomicid', 'strand', 'start', 'end', 'transcript_id'],
            stop_codon => ['genomicid', 'strand', 'start', 'end', 'transcript_id'],
            transcript => ['transcriptid', 'geneid']

        },
        tidcounter => 1,
        genecounter => 1,
        ignore_tag => ['STS'],
        debug => 1
    };
    foreach my $f (keys %{$ref->{column}}) {
        open my $fh, ">$f.tab" or die $!;
        print $fh join("\t", @{$ref->{column}{$f}}), "\n";
        $ref->{$f} = $fh;
    }
    bless ($ref, $class);
    return $ref;
}

sub debug {
    my $self=shift;
    $self->{debug}=1;
}

=head2 contigToTable

contigs from Broad is the key for the GTF files column1

the format is contig_1.409
contig_1. probably means the version.
the number after the period is the contig number.

=cut
sub contigToTable {
    my $self=shift;
    my $file=shift;
    my ($seq, $seqid, $outfile);
    if (!$file) { $file="contigs.fasta"; }
    print STDERR "processing contig file $file ...\n";
    open IN, "<$file" or die $!;
    $outfile=$file;
    $outfile =~ s/.fasta/.tab/;
    open OU, ">$outfile" or die $!;
    $_=<IN>;
    while ($_ && /^>contig_1\.(\d+) /) {
        chomp;
        $seqid=$1;
        $seq="";
        $_=<IN>;
        while ($_ && !/^>/) {
            chomp;
            $seq .= $_;
            $_=<IN>;
        }
        print OU "$seqid\t$seq\n";
    }
    print STDERR "contig file $file converted tabular format\n";
}

=head2 proteinToTablee

    Write seqid title pepseq 
    in a tabular format into $file.tab where $file is the input
    file name given to this function as argument.

=cut
sub proteinToTable {
    my $self=shift;
    my $file=shift;
    my ($seq, $seqid, $title);
    open IN, "<$file" or die $!;
    open OU, ">$file.tab" or die $!;
    $_=<IN>;
    while ($_) {
        chomp;
        if (!/^>/) { die "expecting header line of fasta format\n"; }
        if (/>CC1G_(\d+) \| (.+)/) {
            $seqid=$1;
            $title=$2;
            $title =~ s/ \(translation\) \(\d+ aa\)//;
            $title =~ s/Coprinus cinereus //;
        }
        else { die "unexpted Broad protein id\n"; }
        $seq="";
        $_=<IN>;
        while ($_ && !/^>/) {
            chomp;
            $seq .= $_;
            $_=<IN>;
        }
        $seqid =~ s/^0+//;
        print OU "$seqid\t$title\t$seq\n";
    }
}

=head2 processSequenceFiles
    
    will process all protein and genomic files

    This function is based on the file suffix.
    For NCBI the default is usually 

    genomic: *.gff
    protein: *.fna

    Should have headers as the first line.

    The mit broad institue has only one file
    for each category.

=cut
sub sequenceToTable {
    my $self=shift;
    my $genomic_suffix = shift;
    my $protein_suffix= shift;

    my @files=glob("*.$genomic_suffix");
    my $counter=0;
    my (@arr);
    my ($seq);
    open GNM, ">genomicseq.tab" or die $!;
    print GNM "gi\tgenomicid\ttitle\tsequence\n";
    foreach my $f (@files) {
        if ($self->{debug}) {
            print STDERR "working on file $f ...\n";
        }
        ++$counter;
        open IN, "<$f" or die "Failed to open genomic file $f\n";
# assume one genomic seq per file!
        $_=<IN>; # header
        if (!/^>gi\|(\d+)\|ref\|(.+)\| (.+)/) {
            die "Not proper refseq genomic Fasta format\n";
        }
        @arr=($1,$2,$3);
        $seq="";
        $_=<IN>;
        while ($_ && !/^>/ && !/^\s*$/) {
            chomp;
            $seq .= $_;
            $_=<IN>;
        }
        if ($_ && /^>/) {
            die "more than one genomic sequence per file. Write new code\n";
        }
        print GNM join("\t", @arr), "\t$seq\n";
    }
    print "$counter genomic sequences processed into genomicseq.tab\n";
# protein sequence files, more than one protein per file
    @files = glob("*.$protein_suffix");
    $counter=0;
    open OU, ">proteinseq.tab";
    print OU "gi\tproteinid\ttitle\tequence\n";
    foreach my $f (@files) {
        open IN, "<$f" or die $!;
        $_=<IN>;
        while ($_) {
            chomp;
            if (!/^>/) {
                die "expecting the header line of fasta\n";
            }
            s/\[.+\]//g;
            s/\s$//;
            if (!/^>gi\|(\d+)\|ref\|(.+)\| (.+)/) {
                die "wrong format of protein refseq fasta header\n$_\n";
            }
            my @arr=($1, $2, $3);
            ++$counter;
            $seq="";
            $_=<IN>;
            while ($_ && !/^>/ && !/^\s*$/) {
                chomp;
                $seq .= $_;
                $_=<IN>;
            }
            print OU join("\t", @arr), "\t$seq\n";
        }
    }
    print STDERR "$counter protein sequences processed into proteinseq.tab\n";
}

=head2 processFile($infile)

    This method process one GFF file $infile and write the output
    to filehandles stored in this object. 

    This method will only process GFF files. Separate methods 
    should be called to process the protein and genomic files.

=cut
sub processFile {
    my $self=shift;
    my $file=shift;
    if ($self->{debug}) {
        print STDERR "working on input file $file ...\n";
    }
    open IN, "<$file" or die "$!, Failed to read GTF3 nput file $file\n";

    $_=<IN>;
    my (@arr);
    #my (%header, %t2g);
    my (%t2g);
# %t2g is transcriptid => geneid mapping 
    my ($fh, $ln, $organism, $geneid, $tid, $attr);
    while ($_) {
        chomp;
        if ($_ && (/^\s*$/ || /^#/)) {
            $_=<IN>;
            next;
        }
        @arr=split /\t/;
# clean up junk from Broad
        # supercont_1.198%20of%20Coprinus%20cinereus
        if ($arr[0] =~ /supercont_\d\.(\d+)/) {
            $arr[0] = $1;
        }
        else {
            die "unexpected Broad format\n";
        }
        $attr=makeAttributeHash($arr[8]);
        if (exists $self->{$arr[2]}) { $fh=$self->{$arr[2]}; }
        elsif (grep $arr[2] eq $_,  @{$self->{ignore_tag}}) {
            #print STDERR "Ignoring feature $arr[2]\n";
            $_=<IN>;
            next;
        }
        else {
            die "$arr[2] found in file $file  no associated file handle\n";
        }

        if ($arr[2] eq 'CDS') {
            print $fh $arr[0], "\t", $arr[6], "\t", $arr[7], "\t", $arr[3], "\t", $arr[4];
            if (!exists $attr->{transcript_id}) {
                die "must have transcript_id attribute for CDS feature\n";
            }
            $tid=$attr->{transcript_id};
            delete $attr->{transcript_id};
            print $fh "\t$tid\n"; 
            if (exists $attr->{gene_id}) {
                $geneid=$attr->{gene_id};
                delete $attr->{gene_id};
            }
            else {
                die "CDS must have gene_id attribute\n";
            }
            if (keys %$attr > 0) {
                print STDERR attribute2string($attr), "\n";
                die "attribute other than transcript_id and gene_id\n";
            }
            #$rest=attribute2string($attr);
            transcriptGeneMapping(\%t2g, $tid, $geneid);
        }
        elsif ($arr[2] eq 'start_codon' || $arr[2] eq 'stop_codon'
            || $arr[2] eq 'exon') { 
            if (!exists $attr->{transcript_id}) {
                die "no transcript_id attribute\n";
            }
            $tid=$attr->{transcript_id};
            if (!exists $attr->{gene_id}) {
                die "no gene_id attribute\n";
            }
            $geneid=$attr->{gene_id};
            print $fh $arr[0], "\t", $arr[6], "\t", $arr[3], "\t", $arr[4], "\t$tid\n";
            delete $attr->{transcript_id};
            delete $attr->{gene_id};
            if (keys %$attr > 0) {
                die "attribute not considered", attribute2string($attr), "\n";
            }
            transcriptGeneMapping(\%t2g, $tid, $geneid);
        }
        else {
            die "new feature $arr[2] not considered\n";
        }
        $_=<IN>;
        while ($_ && /^#/) { $_=<IN>; }
    }
# transcript output
    my ($t, $g);
    $fh=$self->{transcript};
    while (($t,$g) = each %t2g) {
        print $fh "$t\t$g\n"
    }
    print "writting to transcript 2 gene mapping file done\n";
}

# helper function
sub transcriptGeneMapping {
    my $h=shift;
    my $tid=shift;
    my $gid=shift;
    if (exists $h->{$tid}) {
        if ($h->{$tid} ne $gid) {
            die "transcriptid => geneid inconsistency from CDS\n";
        }
    }
    else {
        $h->{$tid}=$gid;
    }
}

sub showOrganismInformation {
    my $self=shift;
    my $fh=shift;
    if (!$fh) { $fh = *STDERR; }
    print $fh "Organism_name\tattributes\n";
    my $org=$self->{organism};
    foreach my $o (keys %$org) {
        print $fh $o, "\t", $org->{$o}, "\n";
    }
}

sub outputRest {
    my $attr=shift;
    my $fh=shift;
    my $rest=attribute2string($attr);
    print $fh "\t";
    if ($rest) { print $fh "$rest"; }
    print $fh "\n";
}

# GTF used a different format
# a GTF attribute field
# gene_id "CC1G_00001"; transcript_id "CC1T_00001";
sub makeAttributeHash {
    my $txt=shift; # the 9th column of GFF3
    if (!$txt || $txt =~ /^\s*$/) {
        #warn "empty attribute\n";
        return undef;
    }
    $txt =~ s/;$//;
    while (my ($k,$v) = each %HTMLENCODING) {
        $txt =~ s/$k/$v/g;
    }
    my @arr= split /;\s*/, $txt;
    my ($i, $j);
    my %h;
    for ($i=0; $i<@arr; $i++) {
        if ($arr[$i] =~ /([_a-z]+) "(.+)"/) {
            $h{$1}=$2;
        }
        else {
            die "Wrong attribute format $arr[$i]\n";
        }
    }
    return \%h;
}

=begin comment

sub makeAttributeHash {
    my $txt=shift; # the 9th column of GFF3
    if (!$txt || $txt =~ /^\s*$/) {
        #warn "empty attribute\n";
        return undef;
    }
    my @arr= split /;\s*/, $txt;
    my ($i, $j);
    my %h;
    for ($i=0; $i<@arr; $i++) {
        $j=index($arr[$i], '=');
        if ($j == -1) {
            die "Attribute should have = sign\n";
        }
        if ($j == length($arr[$i])-1) {
            $h{substr($arr[$i], 0, $j)}=1;
        }
        else {
            $h{substr($arr[$i], 0, $j)} = substr($arr[$i], $j+1);
        }
    }
    return \%h;
}

=end comment

=cut

sub attribute2string {
    my $h=shift;
    my @arr;
    my ($k,$v);
    while (($k,$v)= each %$h) {
        push @arr, "$k=$v";
    }
    return join('; ', @arr);
}

# for debug
sub showHash {
    my $h=shift;
    return if (!$h);
    foreach my $k (keys %$h) {
        print STDERR $k, "=>", $h->{$k}, "\n";
    }
}

1;
