#!/usr/bin/perl -w

use strict;

# simply nr sequence id from NCBI use only the primary one
#>gi|295901438|dbj|AB520752.1| Gonium viridistellatum genes for ITS1, 5.8S rRNA
#to 
#>gi295901438
my ($id,$def);
while (<>) {
    chomp;
	if (/^>/) {
		if (/^>([a-z]{2,3})\|(\d+)(.+)/) {
			$id=$1 . $2;
            $def=trimleading($3);
        }
        elsif (/^>([A-Z]{1,5}\d+)(.+)/) {
            $id=$1;
            $def=trimleading($2);
        }
        else {
            die $_, "\nnot considered\n";
        }
        print ">$id $def\n";
	}
    elsif (/^\s*$/) { 
# ignore
    }
    else {
        print $_, "\n";
    }
}

sub trimleading {
    my $str=shift;
    $str =~ s/^\|//;
    $str =~ s/^ //;
    return $str;
}
