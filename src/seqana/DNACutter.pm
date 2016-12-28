#!/usr/bin/perl -w

#this package is for Web application

package DNACutter;
use Exporter;
our @ISA=("Exporter");
our @EXPORT=qw( $HEADER );
use bioseq;
use SVGUtil;
use strict;

=item
	my $cutter = DNACutter->new;
	$cutter->choseEnzymes($choice);
	# here $choice can be "ALL", "GT\d", "\d", "TYPE\d" or \@ARRAY of enzymes
	$cutter->getseq($seqstr);
	# here seqstr can be in fasta, or genbank format, or any string

	Display methods
	$cutter->showbyseq($width, $showruler, $showlower);
	$cutter->showbyenz();
	$cutter->showbyenz_web();
	$cutter->showbypos();
	$cutter->showbypos_web();
	$cutter->showbygraph();

=cut

# use hash to store enzyme information
#                     0            1          2             3
# enzyme_name => [perl_pattern, cut_site, $enzyme_type, $recog_base_cnt ]
our %cutterinfo;
# this vairable should be loaded only once during the 
# http start up
our @field_anno;

our %cuttertab;  # cutter info from 
#ALL ENZYMES (INDIVIDUALLY REFERENCED) WITH ISOSCHIZOMERS

our %suppliertab;

# use global variables faster than into a package
# ?
#my %byname=();
#my %bypos=();
#my @nocut=();
# initialize global variables
if (! %cutterinfo) {
	# when put into the web server it needs absolute path
	read_proto('/usr/local/apache2/perlp/bioinfo/restrict.site');
	read_data('/usr/local/apache2/perlp/bioinfo/restrict.dat');
}
our $HEADER = <<ENDH;
<head><title>Restriction Site Mapping of DNA Sequences</title>
<meta name="keywords" content="restriction enzyme, DNA digestion, nucleic acids, cutter">
<link rel=stylesheet type="text/css" href=/css/contextmenu.css>
<link rel=stylesheet type="text/css" href=/kz.css>
<script src="/jscript/restrict.js"></script>
</head>
ENDH

sub new {
	my $invocant = shift;
	my $class = ref($invocant) || $invocant;
	my $self = { };
	bless($self, $class);
	return $self;
}

# enzyme_name, recognition_string, cut_size, type (1,2,3)
#    0               1                  2          3         
# $enzyme_name, toregex($recogsite), $cutsite, $enzyme_type

#@resite=();
#read_proto($sitefile);
# restriction site information
#  
# enzyme_name, recognition_string, cut_size, type (1,2,3), recog_base_count
#readseq($seqfile);

# put results into two containers
# 1. by enzyme name, 2. by position

# choice can be any digit [1-9], ALL, 6MORE
sub choseEnzymes {
	my $self = shift;
	my $choice = shift;  # ALL or reference to array of enzymes from list
	my %ch = ();
	my ($E, $L);
	if ($choice eq 'ALL') {
		$self->{choice} = \%cutterinfo;
	}
	elsif ($choice =~ /TYPE(\d)/) {
		# only allowed the values are 1,2, and 3
		$L = $1;
		foreach $E (keys %cutterinfo) {
			if ($cutterinfo{$E}->[2] == $L) {
				$ch{$E} = $cutterinfo{$E};
			}
		}
		$self->{choice} = \%ch;
	}
	elsif ($choice =~ /GT(.+)/) {
		# copy 6 or more cutters
		$L = $1;
		foreach $E (keys %cutterinfo) {
			if ($cutterinfo{$E}->[3] > $L) {
				$ch{$E} = $cutterinfo{$E};
			}
		}
		$self->{choice} = \%ch;
	}
	elsif ($choice =~ /^(\d)$/) {
		$L = $1;
		foreach $E (keys %cutterinfo) {
			if ($cutterinfo{$E}->[3] == $L) {
				$ch{$E} = $cutterinfo{$E};
			}
		}
		$self->{choice} = \%ch;
	}
	elsif (ref($choice) eq 'ARRAY') {
		# pick the subset of enzymes from the list
		# by slicing of the hash
		#print "Given an Array!\n";
		foreach $E (@$choice) {
			$ch{$E} = $cutterinfo{$E};
			#print "Enzyme $E ", join('===', @{$cutterinfo{$E}}), "\n";
		}
		$self->{choice} = \%ch;
	}
	else {
		die "$choice is invalid\n";
	}
	#$self->cutseq;
}

# this is the core algorithm to cut the sequence
sub cutseq {
	my $self = shift;
	my $cr = $self->{choice};  # reference to enzymes to use
	my %byname=();
	my %bypos=();
	my @nocut=();
	my $seq = $self->{sequence};
	foreach my $enz (keys %$cr) {
		my $pat = $cr->{$enz}->[0];  # recognistion pattern in perl regex
		my @cutpos=();
		while ($seq =~ /$pat/g) {
			my $position = pos($seq) - length($&);
			push @cutpos, $position;
			if (exists $bypos{$position} ) {
				push @{$bypos{$position}}, $enz;
			}
			else { # does not exist
				$bypos{$position} = [$enz];
			}
		}
		if (@cutpos) {
			$byname{$enz} = [ $cr->{$enz}->[1], @cutpos ];
		}
		else {
			push @nocut, $enz;
		}
	}
	$self->{bypos} = \%bypos;
	$self->{byname} = \%byname;
	$self->{nocut} = \@nocut;
}

sub showbygraph {
	my $self = shift;
	my $W = shift; 
	if (!$W) { $W = 900; } # width of the graph
	my $L = length($self->{sequence});
	my $scale = $W/$L;
	my $posh = $self->{bypos};

	# remap position into a short range
	my %newpos=();
	my $p;
	# no need to sort here
	foreach $p (keys %{$posh}) {
		my $np = (int($scale*$p*10 + 0.5))/10;
		if (exists $newpos{$np}) {
			push @{$newpos{$np}}, @{$posh->{$p}};
		}
		else {
			$newpos{$np} = [ @{$posh->{$p}} ];
		}
	}

	my ($xpos, $graph_buf, $text_buf);
	# keep the last line and position occupied by text here
	# (ypos, xbegin, xend)
	my @prestack=();
	my $lnsp = 13;

	$graph_buf = "";
	$text_buf = "";
	my $maxEnd = 0;  # max end

	foreach $p (sort { $a <=> $b } (keys %newpos)) {
		my $enzs = join(', ', @{$newpos{$p}});
		my $enzslen = length($enzs)*4.8; # num pixels per char

		# first clean the stacks
		# any line further than $lnsp away will be discarded
		# remove far away elements
		while (@prestack && $p - $prestack[0]->[0] > $lnsp) {
				shift @prestack;
		}

		# sort the ranges and find hole enough to hold new segment
		my @segs = ();
		for (my $i=0; $i<@prestack; $i++) {
			push @segs, [ $prestack[$i]->[1], $prestack[$i]->[2] ];
		}

		if (scalar(@prestack) < 1) {
			$xpos = 10;
		}
		else {
			@segs = sort { $a->[0] <=> $b->[0] } @segs;
			if ($segs[0]->[0] > $enzslen + 30) {
				$xpos = 10;
			}
			else {
				$xpos = 0;
				for (my $s=0; $s<$#segs; $s++) {
					if ($segs[$s+1]->[0] - $segs[$s]->[1] > $enzslen + 20) {
						$xpos = $segs[$s]->[1] + 10;
						last;
					}
				}
				if (!$xpos) { $xpos = $segs[$#segs]->[1] + 12; }
			}
		}

		$graph_buf .= "<line x1='0' y1='$p' x2='$xpos' y2='$p' style='stroke:#33ff33;stroke-width:0.5' />\n";
		$text_buf .= "<text x='$xpos' y='". $p . "' style='fill:red;font-size:11;font-family:serif'>$enzs</text>\n";

		push @prestack, [$p, $xpos, $xpos+$enzslen];
		if ($xpos + $enzslen > $maxEnd) { $maxEnd = $xpos+$enzslen; }
	}

	$xpos = 10;
	my $ypos = $maxEnd + 50;
	my $rulery = 80;
	my $jscript = navigatorScript($xpos, $ypos, $rulery);

	my $graph = <<ENDG;
<script><![CDATA[
$jscript
]]> </script>

<text x="20" y="20">Length of sequence: $L</text>
<g transform="translate($xpos, $ypos)" id="main">
<g transform="rotate(-90)">
$graph_buf
<line x1="0" y1="0" x2="0" y2="$W" style="stroke:black;stroke-width:4;" />
$text_buf
</g></g>
ENDG
	print $graph;
	drawNavigator(400, 10);
	print "<g transform='translate($xpos, $rulery)' id='ruler'>\n" . drawRuler($L, 20, $scale) . "</g>"; 
}

sub showbypos {
	my $self = shift;
	my $bypos_r = $self->{bypos};

	foreach my $p (sort {$a <=> $b } keys %$bypos_r ) {
		print $p, "\t", join(', ', @{$bypos_r->{$p}}), "\n";
	}
}
sub showbypos_web {
	my $self = shift;
	my $bypos_r = $self->{bypos};

	print "<tr bgcolor=#ffeeff><th>Position</th><th>Enzymes that Cut</th></tr>\n";
	foreach my $p (sort {$a <=> $b } keys %$bypos_r ) {
		print "<tr><td>$p</td><td>", join(', ', @{$bypos_r->{$p}}), "</td></tr>\n";
	}
}

sub showbyseq {
	my $self = shift;
	my $len = shift;  # width of the line
	my $showruler = shift;
	my $showlower = shift;
	if (!$len) { $len = 70; }
	my $bypos_r = $self->{bypos};
	my $seq = $self->{sequence};
	my $seqlen = length($seq);
	my $linestr = ' ' x ($len+20); # empty line
	#my $ruler = '-' x $len;
	#for (my $i=1; $i<=$len; $i++) {
	#	if ($i % 10 == 0) {
	#		substr($ruler, $i-1, 1) = '+';
	#	}
	#}

	my $p = 0;
	while ($p<$seqlen) {
		my @lines=();
		push @lines, $linestr;

		my $i=0;
		my $begin = $p;
		while ($p<$seqlen && $i<$len ) {
			if (exists $bypos_r->{$p}) {
				# all enzymes at this position
				my @arr = @{$bypos_r->{$p}};
				foreach my $e (@arr) {
					my $wordlen = length($e);
					# look for existing empty lines
					my $found_space = 0; # assume no space
					my $j;
					for ($j = 0; $j<@lines; $j++) {
						my $space = ' ' x ($wordlen+1);
						if (substr($lines[$j], ($i-1), ($wordlen+1)) eq $space ) {
							#substr($lines[$j], $i, $wordlen) = addwebinfo($e, $p);
							substr($lines[$j], $i, $wordlen) = $e;
							$found_space = 1;
							last;
						}
					}
					# add another line
					if (!$found_space) {
						push @lines, $linestr;
						#substr($lines[$j], $i, $wordlen) = addwebinfo($e, $p);
						# cause trouble in this
						substr($lines[$j], $i, $wordlen) = $e;
					}
				}
			}
			$i++;
			$p++;
		}
		# print result
		for ($i=$#lines; $i>=0; $i--) {
			print $lines[$i], "\n";
		}
		my $seqline = substr($seq, $begin, $len);
		#print substr($seq, $begin, $len), "  $p\n";
		print $seqline, "  $p\n";
		if ($showruler) {
			#my $rr = $ruler;
			my $ruler = '-' x $len;
			my $r = 0;
			for (my $i=$begin; $i<$p; $i++) {
				if (($i+1) % 10 == 0) {
					substr($ruler, $r, 1) = '+';
				}
				++$r;
			}
			print "$ruler\n";
		}
		if ($showlower) {
			print complement($seqline), "\n";
		}
		print "\n";
	}
}
# help function
sub addwebinfo {
	my $enz = shift;
	my $info = shift;
	my $id = $enz . "_postip";
	my $webstr = <<ENDW;
<span onmouseover="showtip('$id')" onmouseout="hidetip('$id')" 
	style="cursor:pointer;color:green;">$enz</span> 
<div id='$id' class="contextMenus">$info</div>
ENDW
	return $webstr;
}

sub showbyenz {
	my $self = shift;
	my $byname_r = $self->{byname};

	foreach my $e (sort keys %$byname_r) {
		print $e, "\t", join(', ', @{$byname_r->{$e}}), "\n";
	}
	print "\nEnzymes that did not cut\n";
	print join(', ', @{$self->{nocut}}), "\n";
	print "\n";
}
sub showbyenz_web {
	my $self = shift;
	my $byname_r = $self->{byname};

	print "<tr><th>Enzyme</th><th>Recognition Site</th><th>Location</th></tr>\n";
	foreach my $e (sort keys %$byname_r) {
		my @infoarr = @{$byname_r->{$e}};
		print "<tr><td onclick=\"flipinfo('$e')\" class='flipper'>$e</td><td>", $infoarr[0], "</td><td>";
		print join(', ', @infoarr[1..$#infoarr]), "</td></tr>\n";
		# info for this enzyme
		print "<tr><td colspan=3><div id='$e' class='expandinfo'><table border=1 bgcolor=#ffeeaa>";
		my $info = $cuttertab{$e};
		for (my $i=0; $i<@$info; $i++) {
			if ($i == 5) {
				print "<tr><td><a href='/bioinfo/restrictsup.html' target='sup'>", $field_anno[$i+1]->[0], "</a>";
			}
			else {
				print "<tr><td>", $field_anno[$i+1]->[0];
			}
			print "</td><td>", $info->[$i], "</td></tr>\n";
		}
		print "</table></div></td></tr>\n";
	}
	print "<tr><td colspan=3 bgcolor=#ffeedd>Enzymes that did not cut</td></tr>\n",
		"<tr><td colspan=3 bgcolor=#ffffbb>",
		join(', ', @{$self->{nocut}}), "</td></tr>\n";
}

# object method $obj->readseq($seqfile)
# sequence file in fasta format or just the sequences
sub readseq {
	my $self = shift;
	my $seqfile = shift;
	my $seq="";
	open IN, "<$seqfile" or die "Could not open $seqfile $!\n";
	$_=<IN>;
	chomp;
	if (/^>/) { $self->{title} = $_; }
	else { $seq = $_; }

	while (<IN>) {
		chomp;
		$seq .= $_;
	}
	#print length($seq), "\n";
	close IN;
	$self->{sequence} = $seq;
	if (! exists $self->{choice} ) {
		$self->choseEnzymes('GT5');
	}
	$self->cutseq;
}

sub getseq {
	my $self = shift;
	my $seqstr = shift;
	if ($seqstr =~ /^>(.+)\n/g) {
		$self->{title} = $1;
		my $pos = pos($seqstr);
		$seqstr = substr($seqstr, $pos);
	}
	$seqstr =~ s/[^A-Za-z]//g;
	$self->{sequence} = $seqstr;
	if (! exists $self->{choice} ) {
		$self->choseEnzymes('GT5');
	}
	$self->cutseq;
}

# this is an object method
sub read_proto {
	#my $self = shift;
	my $sitefile = shift;
	my ($enzyme_name, $enzyme_type);

	# read restriction site use the proto format
	open IN, "<$sitefile" or die "Could not open $sitefile $!\n";
	$_ = <IN>;
	while ($_) {
		if (/^#/ | /^\s*$/) { $_ = <IN>; next; }
		if (/TYPE (\d) ENZYMES/) {
			$enzyme_type = $1;
			last;
		}
	}
	# now have read the header for each group of enzyme
	$_=<IN>;
	while ($_) {
		while ($_ && !/TYPE \d ENZYMES/) {
			if (/^#/ | /^\s*$/) { $_=<IN>; next; }
			chomp;
			/(\w+)\s+([-A-Z\^()\d\s\/]+)/;
			$enzyme_name=$1;
			my $cutsite = $2;
			my $recogsite = $cutsite;
			$recogsite =~ s/[\s(\d\/()-]//g;
			$recogsite =~ s/\^//g;
			# we will put recogsite into perl regular expression format
			my $tmp=revcomp($recogsite);
			my $pat = toregex($recogsite);
			if ($tmp ne $recogsite) {
				$tmp = toregex($tmp);
				$pat .= "|$tmp";
			}
			$pat = qr/$pat/i;
			$cutterinfo{$enzyme_name} = [ $pat, $cutsite, $enzyme_type, countRecogsite($recogsite) ];
			$_=<IN>;
		}
		last if (!$_);
		if (/TYPE (\d) ENZYMES/) {
			$enzyme_type = $1;
			$_=<IN>;
		}
	}
	close IN;
}

sub countRecogsite {
	my $cnt = 0;
	my @arr = split //, uc($_[0]);
	for (my $i=0; $i<@arr; $i++) {
		my $b = $arr[$i];
		if ($b eq 'A' || $b eq 'C' || $b eq 'G' || $b eq 'T') {
			$cnt++;
		}
		elsif ($b eq 'R' || $b eq 'Y' || $b eq 'M' || $b eq 'K' || $b eq 'S' || $b eq 'W') {
			$cnt += 0.5;
		}
		elsif ($b eq 'B' || $b eq 'D' || $b eq 'H' || $b eq 'V') {
			$cnt += 1/3;
		}
	}
	return $cnt;
}

# helper function to transform into perl regular expression
sub toregex {
	my $recogsite = shift;
	# we will put recogsite into perl regular expression
	# format
=item
	R = G or A
	Y = C or T
	M = A or C
	K = G or T
	S = G or C
	W = A or T
	B = not A (C or G or T)
	D = not C (A or G or T)
	H = not G (A or C or T)
	V = not T (A or C or G)
=cut

	$recogsite =~ s/N/\./g;
	$recogsite =~ s/R/\[GA\]/g;
	$recogsite =~ s/Y/\[CT\]/g;
	$recogsite =~ s/M/\[AC\]/g;
	$recogsite =~ s/K/\[GT\]/g;
	$recogsite =~ s/S/\[GC\]/g;
	$recogsite =~ s/W/\[AT\]/g;
	$recogsite =~ s/B/\[CGT\]/g;
	$recogsite =~ s/D/\[AGT\]/g;
	$recogsite =~ s/H/\[ACT\]/g;
	$recogsite =~ s/V/\[ACG\]/g;
	#$recogsite = qr/$recogsite/i;
	return $recogsite;
}

# this is an object method
sub read_data {
	#my $self = shift;
	my $sitefile = shift;
	@field_anno=();

	# read restriction site use the proto format
	open IN, "<$sitefile" or die "Could not open $sitefile $!\n";
	$_ = <IN>;
	while (!/^</) { $_ = <IN>; }
	while (/^<([A-Z ]+)>(.*)/ && !/^REBASE codes for commercial sources/) {
		my $header = $1;
		my $anno = $2;
		$_ = <IN>;
		while (!/^</ && !/^REBASE codes for commercial sources/) {
			$anno .= $_;
			$_ = <IN>;
		}
		$anno =~ s/^\s+//;
		push @field_anno, [$header, $anno];
		#print "$header => $anno\n";
	}
# read rebase comercial code
	$_=<IN>;
	while (/^\s*$/) { $_ = <IN>; }
	while (!/^<1>/ && !/^\s*$/) {
		chomp;
		/^\s+([A-Z])\s+(.+)/;
		$suppliertab{$1}=$2;
		#print "Supplier: $1 => $2\n";
		$_=<IN>;
	}
	while (!/^<1>/) { $_ = <IN>; }
=item
there is no need to compile the perl pattern from
this file, because this file lists multiple recognition
site for each enzyme, just use this file as informational
=cut

	while ($_ && /^<1>/) {
		chomp;
		my $enz = substr($_, 3);
		my @arr = ();
		for (my $i=2; $i<9; $i++) {
			$_ = <IN>;
			chomp;
			push @arr, substr($_, 3);
		}
		$cuttertab{$enz} = [  @arr ];

=item
		$recogsite = $arr[1];
		$recogsite =~ s/[\s(\d\/()-]//g;
		$recogsite =~ s/\^//g;
		# we will put recogsite into perl regular expression format
		#print "$enz Recogsite: $recogsite\n";
		# enzyme does not have recognition site, useless to this program
		if ($recogsite ne '?') {
			my $tmp=revcomp($recogsite);
			my $pat = toregex($recogsite);
			if ($tmp ne $recogsite) {
				$tmp = toregex($tmp);
				$pat .= "|$tmp";
			}
			print $pat, "\n";
			$pat = qr/$pat/i;
			$cuttertab{$enz} = [ $pat, @arr, countRecogsite($recogsite) ];
		}
=cut

		$_ = <IN>;
		while ($_ && !/^<1>/) { $_ = <IN>; }
	}

	close IN;
}

sub dumpCuttertab {
	open OU, ">cutter.tab";
	foreach my $e (keys %cuttertab) {
		print OU "$e\t", join("\t", @{$cuttertab{$e}}), "\n";
	}
	close OU;
}

sub dumpSuppliertab {
	open OU, ">supplier.tab";
	foreach my $s (keys %suppliertab) {
		print OU "$s\t", $suppliertab{$s}, "\n";
	}
	close OU;
}

sub getEnzymeList {
	my $self=shift;
	return sort(keys %cutterinfo);
}

# print form
sub printForm {
	my $self=shift;
#<a href=/index.html><img height=36 src=/icons/logo.gif border=0></a>
#<span style="font-family:cursive;font-size:400%;font-weight:bold;" >
#	DNA Cutter</span>
# these should be done on the driver program

print <<ENDF;
<form method="POST" action="/perlp/bioinfo/cut" target="resmap" style="display:inline">
<table bgcolor=#bbeeff>
	<tr><td colspan=3>
		<span onmouseover="showtip('seqfmtip')" onmouseout="hidetip('seqfmtip')" style="cursor:pointer">Cut and Paste Sequence into the Text Box</span>
		<div id="seqfmtip" class="contextMenus">
			Input sequence can be in most formats:<br>
<pre>
Fasta:
>seqname title
ACGTAAAATTTTTTT
GenBank
 1 cgccgcggga ggcggacgag atgcgagcgc ...
61 tgctggcgct gggggcgctg gcgggcgttg ...
</pre>
		</div>
		</td></tr>
	<tr><td colspan=3>
<textarea name="seq" cols=90 rows=10>
</textarea></td></tr>
<tr><td>Enzymes</td>
	<td colspan=2>
	<input type=radio name=choice value="ALL"> All
	<input type=radio name=choice value="TYPE1"> Type I
	<input type=radio name=choice value="TYPE2"> Type II
	<input type=radio name=choice value="TYPE3"> Type III
	<input type=radio name=choice value="bysite" checked> 
		<span id="recog" onmouseover="showtip('recogtip')" onmouseout="hidetip('recogtip')" style="cursor:pointer">Recogsite</span>
		<div id="recogtip" class="contextMenus">
			Computed as the sum for each Base.<br>
			If all bases in the recognition site are<br>
			unambiguous then it is the length of the<br>
			recognition site.<br>
			When Counting bases:<br>
			A, C, G, or T count 1. <br>
			Ambiguous positions:<br>
			R = G or A; Y = C or T<br>
			M = A or C; K = G or T<br>
			S = G or C; W = A or T<br>
			count 1/2 <br>
			B = C or G or T; D = A or G or T<br>
			H = A or C or T; V = A or C or G<br>
			count 1/3 and N count 0
		</div>
		&gt;  <input type=text name=recogsize size=3 value="5.9"></td></tr>
<tr valign=top><td></td>
	<td valign=top colspan=2>
	<input type=radio name=choice value="custom"> Pick From List 
	<select name=enzlist multiple="multiple" size=4 title="Type a letter to jump to Enzymes starting with that letter. Hold-on Ctrl button to select multiple enzymes">
ENDF

	foreach my $e ( $self->getEnzymeList ) {
		print "<option>$e</option>\n";
	}

print <<ENDF;
	</select>Press Ctrl Key for Multiple Choice</td></tr>
<tr><td>Display</td>
	<td><input type=text name=width value=80 size=3>nt per line</td>
	<td>Ruler: <input type=radio name=showruler value=1 checked> yes
	          <input type=radio name=showruler value=0> no
	</td></tr>
<tr><td colspan=2>
		Lower Strand:
				<input type=radio name=showlower value=1> yes
	          <input type=radio name=showlower value=0 checked> no
	</td>
</tr>
</table>
Draw <input type=submit name=output value=Text onclick="focus_map()"> 
	Graph Width <input name=graph_width value=900 size=4> <input type=submit name=output value=Graph onclick="focus_map()">
</form>
ENDF
}

sub printHeader {
	print <<ENDH;
<head><title>Restriction Site Mapping of DNA Sequences</title>
<link rel=stylesheet type="text/css" href=/css/contextmenu.css>
<link rel=stylesheet type="text/css" href=/kz.css>
<script src="/jscript/restrict.js"></script>
</head>
ENDH

}

1;
