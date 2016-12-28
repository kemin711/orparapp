package GFF3Converter;

#use Exporter;
use strict;

#our @ISA=qw(Exporter);
#our @EXPORT=qw(readfile);


=head1 NAME

GFF3Converter - convert GFF3 into Various formats

=head1 DESCRIPTION

I will use this perl module to convert GFF3 formats into mainly
relational database schema which in turn could be easily changed
to other formats.

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
        column => { source => ['genomicid', 'strand', 'start', 'end', 'mol_type', 'chromosome', 'organism'],
                    CDS => ['genomicid', 'strand', 'phase', 'start', 'end', 'id', 'exon_number', 'protein_id'],
                    exon => ['genomicid', 'strand', 'start', 'end', 'exon_number', 'transcript_id'],
                    intron => ['genomicid', 'strand', 'start', 'end', 'id', 'number', 'parent', 'other_attributes'],
                    gene => ['genomicid', 'strand', 'start', 'end', 'id_or_locus_tag', 'attribute'],
                    repeat_region => ['genomicid', 'strand', 'start', 'end', 'attribute'],
                    misc_feature => ['genomicid', 'strand', 'start', 'end', 'attribute'],
                    start_codon => ['genomicid', 'strand', 'start', 'end', 'protein_id'],
                    stop_codon => ['genomicid', 'strand', 'start', 'end', 'protein_id'],
                    protein => ['proteinid', 'geneid', 'product', 'attribute'],
                    transcript => ['transcriptid', 'RNA_type', 'geneid', 'attribute']

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

=head2 processSequenceFiles
    
    will process all protein and genomic files

    This function is based on the file suffix.
    For NCBI the default is usually 

    genomic: *.gff
    protein: *.fna

    Should have headers as the first line.

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
    print OU "\gi\proteinid\ttitle\sequence\n";
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
    open IN, "<$file" or die "$!, Failed to read GFF3 input file $file\n";

    $_=<IN>;
    my (@row, @arr);
    my (%header, %protein, %transcript);
    my ($ln, $organism);
# read header part, so far not useful
    while ($_ && /^#/) {
        chomp;
        $ln = substr($_, 2);
        $ln =~ /(.+?) (.+)/;
        $header{$1}=$2;
        $_=<IN>;
    }
    while ($_) {
        while ($_ && !/^\s*$/ && !/^#/) {
            chomp;
            my $attr;
            @arr=split /\t/;
            if ($arr[8]) {
                while (my ($k,$v) = each %HTMLENCODING) {
                    $arr[8] =~ s/$k/$v/g;
                }
                $attr=makeAttributeHash($arr[8]);
            }
            my $fh;
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
# CDS genomic strand phase start end id exon_number parent
# other_attribute(k=v;)
# NCBI's NC_ have different format from NW_
                # reals sucks!!!
                print $fh $arr[0], "\t", $arr[6], "\t", $arr[7], "\t", $arr[3], "\t", $arr[4];
                if (exists $attr->{ID}) {
# some CDS has ID attribute some don't !!!
# those from NC has ID attribute
                    print $fh "\t", $attr->{ID}; 
                    delete $attr->{ID};
                }
                else {
                    print $fh "\t\\N";
                }
                print $fh "\t", $attr->{exon_number};
                delete $attr->{exon_number};
                my ($pid, $parent, $rest);
                if (exists $attr->{protein_id}) { # required for CDS
                    $pid=$attr->{protein_id};
                    delete $attr->{protein_id};
                }
                else { 
                    die "no proteinid for CDS\n";
                }
                if (exists $attr->{Parent}) {
                    $parent=$attr->{Parent};
                    delete $attr->{Parent};
                }
                elsif (exists $attr->{locus_tag}) {
                    $parent=$attr->{locus_tag};
                    delete $attr->{locus_tag};
                }
                $rest=attribute2string($attr);
                if ($pid) {
                    if (exists $protein{$pid}) {
                        if ($protein{$pid}->[0] ne $parent) {
                            die "different parent (gene) for the same protein\n";
                        }
                        elsif (attribute2string($protein{$pid}->[1]) ne $rest) {
                            die "differnt attribute for the same protein\n";
                        }
                    }
                    else {
                        $protein{$pid}=[$parent, $attr];
                    }
                    print $fh "\t$pid\n";
                }
            }
            elsif ($arr[2] eq 'start_codon' || $arr[2] eq 'stop_codon') { 
                print $fh $arr[0], "\t", $arr[6], "\t", $arr[3], "\t", $arr[4];
                if (exists $attr->{protein_id}) {
                    print $fh "\t", $attr->{protein_id}, "\n";
                }
                else {
                    die "start_codon and stop_codon should have protein_id attribute\n";
                }
            }
            else {
                print $fh $arr[0], "\t", $arr[6], "\t", $arr[3], "\t", $arr[4];
                my ($RNA_type, $tid, $geneid, $rest);
                if ($arr[2] eq 'exon') {
                    print $fh "\t", $attr->{exon_number}, "\t";
                    $RNA_type=$attr->{gbkey};
                    delete $attr->{exon_number};
                    delete $attr->{gbkey};
                    if (exists $attr->{insd_transcript_id}) {
                        $tid=$attr->{insd_transcript_id};
                        delete $attr->{insd_transcript_id};
                    }
                    else {
                        # only mRNA has transcript_id
                        $tid=$self->{tidcounter}++;
                    }
                    print $fh "$tid\n";

                    if (exists $attr->{locus_tag}) {
                        $geneid= $attr->{locus_tag};
                        delete $attr->{locus_tag};
                    }
                    else {
                        $geneid=$self->{genecounter}++;
                    }
                    #print $fh "\t$geneid";
                    if (exists $transcript{$tid}) {
                        if ($transcript{$tid}->[0] ne $RNA_type
                            || $transcript{$tid}->[1] ne $geneid
                            || attribute2string($transcript{$tid}->[2]) ne attribute2string($attr)) {
                            die "same transcript but different attributes\n";
                        }
                    }
                    else {
                        $transcript{$tid}=[$RNA_type, $geneid, $attr];
                    }
                }
                elsif ($arr[2] eq 'gene') {
                    if (exists $attr->{ID}) {
                        print $fh "\t", $attr->{ID};
                        delete $attr->{ID};
                    }
                    elsif (exists $attr->{locus_tag}) {
                        print $fh "\t", $attr->{locus_tag};
                        delete $attr->{locus_tag};
                    }
                    outputRest($attr, $fh);
                }
                elsif ($arr[2] eq 'intron') {
# there is no intron feature for NW_*.gff
                    print $fh "\t", $attr->{ID}, "\t", $attr->{number}, "\t", $attr->{Parent};
                    delete $attr->{ID};
                    delete $attr->{number};
                    delete $attr->{Parent};
                    outputRest($attr, $fh);
                }
                elsif ($arr[2] eq 'source') {
                    print $fh "\t", $attr->{mol_type};
                    delete $attr->{mol_type};
                    if (exists $attr->{chromosome}) {
                        print $fh "\t", $attr->{chromosome};
                        delete $attr->{chromosome};
                    }
                    else {
                        print $fh "\t\\N";
                    }
                    $organism=$attr->{organism};
                    delete $attr->{organism};
                    print $fh "\t$organism\n";

                    if (exists $self->{organism}{$organism}) {
                        if ($self->{organism}{$organism} ne attribute2string($attr)) {
                            die "same organism different attributes!\n";
                        }
                    }
                    else {
                        $self->{organism}{$organism} = attribute2string($attr);
                    }
                }
                elsif ($arr[2] eq 'repeat_region' || $arr[2] eq 'misc_feature') {
# repeat_region usually don't have attribute
                    #print STDERR "repeat_region seen $_\n";
# so far there are only two repeat regions
                    print $fh $arr[0], "\t", $arr[6], "\t", $arr[3], "\t", $arr[4], "\n";
                    if ($attr) {
                        print $fh "\t", attribute2string($attr), "\n";
                        #die "repeat_region got some attribute, need new code\n";
                    }
                    else {
                        print $fh "\t\\N\n";
                    }
                }
                else {
                    die "unseen GFF3 type $arr[2] in file $file\n";
                }
            }
            $_=<IN>;
            while ($_ && /^#/) { $_=<IN>; }
        }
    }
# print out protein
    #print STDERR scalar(keys %protein), " proteins\n";
    my $fh=$self->{protein};
    my ($p, $f);
    while (($p, $f) = each %protein) {
        print $fh "$p\t", $f->[0], "\t"; 
        if (exists $f->[1]->{product}) {
            print $fh $f->[1]->{product}; 
            delete $f->[1]->{product};
        }
        else {
            print $fh "hypothetical protein";
        }
        print $fh "\t", attribute2string($f->[1]), "\n";
    }
# transcript output
    $fh=$self->{transcript};
    while (($p,$f) = each %transcript) {
        print $fh "$p\t", $f->[0], "\t", $f->[1], "\t", attribute2string($f->[2]), "\n";
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
