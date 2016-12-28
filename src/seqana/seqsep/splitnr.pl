### perl ###
##!/usr/bin/perl -w
##!/usr/local/bin/perl -w 

# the nr database is too large, this program will 
# split it into the major division according to the taxonomy database
# when the file size is > 2 GB then this program fails to 
# finish. Currently the nr db is 2.7GB, then the entreis
# after the 2.0 GB are not written to the output file.

use JGIDB;
use DBI;
#use integer;
use Bioseq;
use SeqIterator;
use Taxonomy;
use integer;

my %config=(taxroot => "/house/groupdirs/gat/projects/kzhou/taxonomy",
    nrfile => "/house/blast_db/nrJuly2010/nr",
    action => "nofungi", nsperfile => 4400000);
my %node; # taxnodes
#our (%gi2tax, %taxname2id);
my ($host, $fungiTaxid);
my $a=0;
while ($ARGV[$a]) {
    if ($ARGV[$a] eq '-h') { $host=$ARGV[++$a]; }
    else {
        $config{nrfile}=$ARGV[$a];
    }
    ++$a;
}

if ($config{action} eq 'nofungi') { 
    loadTaxonomy($config{taxroot}, 0, 0);
    nofungi(); 
}
else {
    dividenr();
}

###### sub routines #########

=head2 dividenr

    This is an old implementation with hard coded nr files
    Needs to be rewritten to be more flexible. This method
    also depends on mysql database. Future versions should
    only rely on data files dump from NCBI.  
=cut
sub dividenr {
    my $database='taxonomy';
    my $host='genome-db';
    my $nrfile='/home/blast_db/public_db/nr';
    my $numperfile=3500000; # about 1.8 GB
    my $db=JGIDB->new;
    $db->connect(host=>$host, database=>$database);

    my $sqlstr="select id, name from tax_division";
    my $sth=$db->run($sqlstr);
    my $row;
#my $counter=0;
    while ($row=$sth->fetchrow_arrayref) {
        $name=$row->[1];
        $name =~ s/ /_/g;
        #$div[$row->[0]] = $name;
        my $fh;
        open $fh, ">$name.fas";
        $fhs[$row->[0]]=$fh;
    }
#foreach $d (@div) {
#    print $d, "\n";
#}
    print "The query takes a long time!\n";
    $sqlstr="select p.gi, n.division from prtgi2taxid p join tax_nodes n on p.taxid=n.id";
    $sth=$db->run($sqlstr);
#$counter=0;
    while ($row=$sth->fetchrow_arrayref) {
        #++$counter;
        #if ($counter % 10000 == 0) {
        #    print "at $counter ", $row->[0], "\n";
        #}
        $gi2div{$row->[0]}=$row->[1];
        #my $fh=$fhs[$row->[1]];
        #print $fh $row->[0], "\n";
    }
    print scalar(keys %gi2div), " gi loaded\n";
    open LOG, ">unmapped.fas";

# now read nr file
    my $seqreader=SeqIterator->new($nrfile, 1);
    my $seq=$seqreader->next;
    my ($fh,$div, $i);
    while ($seq) {
        $found=0;
        @ids = $seqreader->getAllIds;
        $id=shift @ids;
        while ($id) {
            if (exists $gi2div{$id}) {
                $found=1;
                last;
            }
            $id=shift @ids;
        }
        if ($found) {
            $div=$gi2div{$id};
            $fh=$fhs[$div];
            print $fh ">$id";
            if (exists $seqreader->{primaryheader}) {
                print $fh " ", $seqreader->{primaryheader};
            }
            print $fh "\n";
            printSeq($seq, $fh);
        }
        else {
            #die "$id not found in gi2div\n";
            print LOG $seqreader->header, "\n";
            $badcnt++;
        }
        $seq=$seqreader->next;
    }
    print $badcnt, " unmapped gi\n";
}

=head2 nofungi

    produce a file of no fungi from nr protein file
    print STDERR "Fungi taxid is

>gi|299754274|ref|XP_001839908.2| hypothetical protein CC1G_06098 [Coprinopsis cinerea okayama7#130] 
gi|298410679|gb|EAU81887.2| hypothetical protein CC1G_06098 [Coprinopsis cinerea okayama7#130]

>299754274 hypothetical protein CC1G_06098 [Coprinopsis cinerea okayama7#130]
This one get into nrnofungi file

=cut
sub nofungi {
    open OU, ">nrNoFungi0.fas" or die $!;

    $fungiTaxid=$taxname2id{Fungi}->[0];
    print STDERR "Fungi taxid is $fungiTaxid\n";
    if ($fungiTaxid != 4751) { # I hard coded this number
        die "fungi taxid has changed, need to update code!\n";
    }
# now read nr file
    my $seqreader=SeqIterator->new($config{nrfile}, 0);
    my $seq=$seqreader->next;
    my ($count, $sp);
    $count=1;
    while ($seq) {
        my $gid = $seqreader->getId;
        #print STDERR "working on $gid\n";
        $sp=$seqreader->getSpecies;
        if (!$sp) { # should not write into output stream
            #warn "header has no species info, discarding\n";
        }
        elsif (nameIsFungi($sp) == 0) {
            if ($count % 200000 == 0) {
                print STDERR "$count sequences processed at $gid ...\n";
            }
            if ($count % $config{nsperfile} == 0) {
                my $vol = $count/$config{nsperfile};
                open OU, ">nrNoFungi$vol.fas" or die $!;
                print STDERR "starting a new file: $vol\n";
            }
            my $title=$seqreader->title;
            if ($title) {
                print OU ">$gid ", $seqreader->title, "\n";
            }
            else {
                print OU ">$gid\n";
            }
            printSeq($seq, *OU);
            ++$count;
            #print STDERR "non fungal sequence $gid\n";
        }
        else {
            #print STDERR "fungal or unknown seq ignored\n";
        }
        $seq=$seqreader->next;
    }
    print STDERR "$count non-fungal sequences written to nrNoFungi*.fas\n";
}

sub nameIsFungi {
    my $n=shift;
    my $tid;
    my $found=0;
    my $shortn;
    if (exists $taxname2id{$n}) {
        $tid=$taxname2id{$n}->[0];
        $found=1;
    }
    if (!$found && ($n =~ / chloroplast$/ || $n =~ / endosymbiont /
                    || $n =~ / x /)
    ) {
        return -1;
    }
    if (!$found && exists $taxname2id{lc($n)}) {
        $tid=$taxname2id{lc($n)}->[0];
        $found=1;
    }

    if (!$found) {
        $shortn=fixName($n);
        if (exists $taxname2id{$shortn}) {
            $tid=$taxname2id{$shortn}->[0];
            $found=1;
        }

        if (!$found && exists $taxname2id{lc($shortn)}) {
            $tid=$taxname2id{lc($shortn)}->[0];
            $found=1;
        }
        else {
            if ($shortn =~ /virus/ || $shortn =~ / sp\./ 
                || $shortn =~ /Enterobacteria phage/
                || $shortn =~ /unidentified/i
                || $shortn =~ /plasmid/i
                || $shortn =~ / expression/i
                || $shortn =~ / vector/i
                ) {
# suppress printing
            }
            elsif ($n ne $shortn) {
                warn "fixed name $shortn not found in hash\n";
            }
            else {
                warn "organism name $n not found in taxname2id hash!\n";
            }
            return -1;
        }
    }
# do look up
    if (rootIs($tid, 4751) == 1) {
        #print "$n is fungi\n";
        return 1;
    }
    else {
        #print "$n is not fungi\n";
        return 0;
    }
}

sub fixName {
    my $n=shift;
    if ($n =~ s/^Synthetic/synthetic/ ||
        $n =~ s/ Virus/ virus/
        || $n =~ s/Cloning vector/cloning vector/
        ) {
    }
    elsif ($n =~ /^[A-Z][a-z]+ sp$/) {
        $n .= '.';
    }
    elsif ($n =~ /virus\b/) {
        $n=upcaseVirusName($n);
    }
    elsif ($n =~ /^([A-Z][a-z]+ [a-z]+) strain/
            || $n =~ /(^[A-Z][a-z]+ [a-z]+) str\. /) {
        $n=$1;
    }
    elsif ($n =~ /^([A-Z][a-z]+ [a-z]+) /) {
        $n=$1;
    }
    elsif ($n =~ /^([A-Z][a-z]+ sp\.) .+/) {
        $n=$1;
    }
    return $n;
}

sub upcaseVirusName {
    my $n=shift;
    if ($n =~ /^([a-z])(.+ virus)$/ 
        || $n =~ /^([a-z])(.+ [a-z]+virus)$/
        || $n =~ /^([a-z])(.+ [a-z]+virus \d+)$/
        || $n =~ /^([a-z])(.+ [a-z]+virus [A-Z])$/
        || $n =~ /^([a-z])(.+ [a-z]+virus [A-Z0-9]+)$/
        || $n =~ /^([a-z])(.+ [a-z]+virus .+)$/
        || $n =~ /^([a-z])([a-z]+virus .+)$/
    ) {
        $n=uc($1) . $2;
    }
    elsif ($n =~ /virus\b/) {
        $n =~ s/^([a-z])/\u$1/;
    }
    return $n;
}

=head2 isFungi

    return 0 if not, 1 yes, -1 not found in database

=cut
sub isFungi {
    my $name=shift;
    my $tid;
    $tid = $taxname2id{$name};
    if (!$tid && ($name =~ s/ strain .+// || $name =~ s/ strain$//
            || $name =~ s/ str\. .+//)) {
        $name =~ s/\s+$//;
        $tid = $taxname2id{$name};
    }
    if ($tid) {
        my $parentid=$tid;
        while (exists $node{$parentid}) {
            if ($fungiTaxid == $parentid) {
                return 1;
            }
            $parentid=$node{$parentid}->[0];
        }
        return 0;
    }
    else {
        #warn "species name: $name not found in tax name2id table\n";
        return -1
    }
}

sub giFromFungi {
    my $gi=shift;
    if (! exists $gi2tax{$gi}) {
        die "$gi not found from gi2tax\n";
    }
    my $taxid=$gi2tax{$gi};
    if ($taxid == $fungiTaxid) { return 1; }
    my $parentid=$node{$taxid}->[0];
    while (exists $node{$parentid}) {
        if ($fungiTaxid == $parentid) {
            return 1;
        }
        $parentid=$node{$taxid}->[0];
    }
    return 0;
}


