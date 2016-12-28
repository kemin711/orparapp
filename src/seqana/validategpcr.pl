#!/usr/bin/perl -w

=for purpose
Valicate new potential gpcr candidates by looking at the title and the 
keywords.  Just my first try.

Iteration object-by-object is very slow.  I have done this way:
1. Dump the keyset of targets
2. Inport into ACEDB
3. Select active keyset, execute aql select :1, :1->title, :1->keyword from @active
4. Dump the result as text.
5. Iterate through ever row in this program.

This is based on title and keyword.
=cut

use dbinfo;
use Pg;
use Ace;

@title_no = ("SH3 domain.binding.protein", 
             "transporting.ATPase",
				 "Leucine.rich.+glycoprotein",
				 "amino.acid.transporter",
				 "ATP.binding.cassette transporter",
				 "Acrosin",
				 "Insulin.like growth factor.binding.protein",
				 "Alu subfamily.+sequence contamination warning entry",
				 "Apolipoprotein", "Arginine metabolism regulation protein",
				 "ATP synthase . chain", "Tyrosine-protein kinase",
				 "actinin", "cytochrome P450",
				 "NADH.ubiquinone oxidoreductase",
				 "Phosphatidylinositol .+kinase", "Toll.like receptor",
				 "Toll protein", "Adenylate cyclase", 
				 "protein disulfide isomerase",
				 " transporter", "Protocadherin", "cotransporter",
				 "Nucleolar phosphoprotein",
				 "agglutinin attachment subunit",
				 "Agrin", "Asporin", "Adenomatous polyposis coli protein",
				 "Aminotriazole resistance protein",
				 "Cadherin.+ precursor",
				 "Brain-cadherin", "Neural-cadherin",
				 "cell adhesion molecule", "Carboxypeptidase",
				 "B.cell receptor","Complement factor",
				 "Chondroadherin", "Chaoptin", "Clathrin",
				 "Complement C\d", "Complement component C\d",
				 "Cartilage oligomeric matrix protein",
				 "Connectin", "Sulfate permease", "Cylicin",
				 "Dentin sialophosphoprotein", 
				 "epidermal growth factor precursor",
				 "Fibulin", "Fibrillin", "Uracil permease",
				 "low.density lipoprotein", "Glucan synthase",
				 "membrane ATPase", "protein disulfide.isomerase",
				 "drug resistance protein", "Heat shock protein",
				 "Chitin synthase", "acetylglucosaminyltransferase",
				 "Mucin", "RNA binding protein", "NADPH thyroid oxidase",
				 "Latrophilin", "Complement receptor",
				 "Polycystic kidney disease", "Metalloprotease disintegrin",
				 "Transcription factor", "Fibropellin", "Fibromodulin",
				 "Keratocan", "Sodium.calcium exchanger",
				 "Protein kinase C-binding protein", "Noelin",
				 "Nidogen", "Osteomodulin", "Opticin",
				 "Biglycan", "Decorin", "Prolargin", "Properdin",
				 "Kit ligand", "Slit protein", "Semaphorin",
				 "Spore coat protein", "Synaptoporin", "Tenascin",
				 "Thoroglobulin", "disulfide isomerase", "Thrombospondin",
				 "NADPH oxidase", "Helicase", "Butyrophilin",
				 "Secreted frizzled-related protein",
				 "Secreted apoptosis related protein",
				 "Thioredoxin related protein", "Tight junction protein",
				 "Thrombospondin-related adhesive protein",
				 "Cytochrome b", "Death receptor Fas", 
				 "Epidermal growth factor"
				);
@title_yes = ("G.protein.coupled.+receptor", 
      "Olfactory receptor", 
		"opsin", 
		"Histamine H\d receptor", "Histamine H. receptor",
		"Histamine receptor",
		"Dopamine receptor D.", "Dopamine .. receptor",
		"Prostaglandin .. receptor", 
		"Chemokine receptor 2", "CC chemokine receptor",
		"Muscarinic acetylcholine receptor",
		"Muscarinic receptor",
		"Calcitonin receptor", "taste receptor",
		"Seven transmembrane helix receptor",
		"seven transmembrane protein",
		"Mu opioid receptor",
		"hydroxytryptamine . receptor",
		"vasopressin receptor"
		);

@keyword_no = ( "Nuclear protein", "DNA.binding",  "Helicase",
		"ATP.binding", "Hydrolase", "Golgi stack", "Protein transport",
		"Actin.binding", "Ionic channel", "Transcription regulation",
		"Mitochondrion", "Kinase", "Blood coagulation", 
		"Oxidoreductase", "Peroxisome", "Purine metabolism", 
		"Transferase", "Plasma", "Complement pathway",
		"cGMP synthesis", "Cytoskeleton", "Structural protein",
		"Myosin", "Lyase", "Sugar transport", "Sodium transport",
      "Lectin", "Silk"
       );

#$source='file'; # default behavior
$delimiter='";"';
$i=0;
if (@ARGV) {
	while ($ARGV[$i]) {
		if ($ARGV[$i] eq '-f') { $file = $ARGV[++$i]; }
		$i++;
	}
}

#@unsure = ();  # Not known whether it is gpcr or not
$swdb = get_acedb("swprt");

if ($file) { fromFile($file); }
else { fromPGDB(); }
#judge_by_motifdb(); # done inside from file method

############################ subroutine #######################
=begin Not_used
sub judge_by_motifdb {
	# this subroutine will judge every entry in @unsure based on
	#InterPro Pfam  PRINTS  SMART PROSITE
	print scalar(@unsure), " negative proteins\n";
	my $swdb = get_acedb("swprt");
	my $tmp;

	my $score;

	foreach $p (@unsure) {
		$score = gpcr_dbmotif($p);
		if ($score) { print "Is GPCR $score\n"; }
		else { print "Not GPCR\n"; }
		#print $p, " is not a GPCR\n\n";
	}
}
=end
=cut

# use the global ACEDB $swdb as source of input
# return the score (number of gpcr_motifs found in all motif db)
sub gpcr_dbmotif {
	my $p = shift;   # input protein identifier
	my $prt= $swdb->fetch(Protein=>$p) or die "Cannot get $p from $swdb\n";
	my $obj = $prt->at("DB_info.DB_XREF");
	if (!$obj) { # no DB_XREF, do nothing
		return 0;  # no motif db entry
		#next;
	}
	$obj = $obj->right;

	#print "\nJudging protein $p\n";
	my $score=0;
	while ($obj) {  # iterate through all referenced dbs
		#print "\nDBXREF: ", $obj, "\n";
		if ($obj eq 'InterPro') {
			$score += visitXREF($obj, "GPCR");
		}
		elsif ($obj eq Pfam) {
			$score += visitXREF($obj, "7tm");
		}
		elsif ($obj eq 'PRINTS') {
			$score += visitXREF($obj, "GPCR");
		}
		elsif ($obj eq 'SMART') {
			$score += visitXREF($obj, "HormR");
		}
		elsif ($obj eq 'PROSITE') {
			$score += visitXREF($obj, "G_PROTEIN_RECEP");
		}
		$obj = $obj->down;
	}
	return $score;
}

sub visitXREF {
	# return true or false (1 or 0)
	my $obj = shift;
	my $pattern = shift;
	#print "\tMotif DB", $obj, "\n";
   $obj = $obj->right;              # traverse the motifs
	my $mtf;

   while ($obj) {
		$mtf = $obj->right->name;
		#print  "\t\t", $mtf, "\n";
		if ($mtf =~ /^$pattern/) {
			return 1;
			#print "\t\t\t is a GPCR ************\n";
		}
		$obj = $obj->down;
	}
	return 0; # false, not GPCR
}

sub fromFile {
	open TOU, ">gpcr.true";
	open FOU, ">gpcr.false";
	open UOU, ">gpcr.unknown";  # hypothetical protein
	open USO, ">gpcr.unsure";  # not sure

	my $inf = shift;
	open IN, "<$inf" or die "Cannot open $inf\n";
	my $ln = <IN>;
	my @arr;
	while ($ln) {
		@arr = toarr($ln);
		$prt = $arr[0];
		if ($arr[1]) { $title = $arr[1]; }
		else { 
			#print "$prt has no title!\n";
			die "$prt has no title!\n";
			$ln = <IN>; 
			next;
		}  # nothing to do

		#print "$prt\n";
		$is_gpcr = title_type($title);
		@keywords = ();                 # initialize to nothing
		while ($prt eq $arr[0] && $title eq $arr[1]) {  # one protein sequence
			if ($arr[2]) { push @keywords, $arr[2]; }
			$ln = <IN>;
			if (!$ln) { last; }
			@arr = toarr($ln);
		}
		if ($is_gpcr > 1) {
			$is_gpcr=keyword_type(\@keywords);
			if ($is_gpcr > 1) {
				$gpcr_motif_cnt = gpcr_dbmotif($prt);
				if ($gpcr_motif_cnt) { $is_gpcr = 1; }
				else { $is_gpcr = 3; }
			}
		}
		output();
		if ($is_gpcr == 0) { print $prt, "\tf\n";  }
		elsif ($is_gpcr == 1) { print $prt, "\tt\n";  }
	}

	close TOU;
	close FOU;
	close UOU;
	close USO;
}

sub show {
	my $fh = shift;
	print $fh $prt, "\n", $title, "\n";
	if (@keywords) { print $fh join('; ', @keywords); }
	print $fh "\n\n";
}

# for debug information
sub output {
	if ($is_gpcr == 1) {
		show(\*TOU);
	}
	elsif ($is_gpcr == 0) {
		#push @unsure, $prt;
		show(\*FOU);
	}
	elsif ($is_gpcr == 2) {
		#push @unsure, $prt;
		show(\*UOU);
	}
	else { # 3 unsure cannot make decision
		#push @unsure, $prt;
		show(\*USO);
	}
}

# 0 false, 1 true, 2 unknown, 3 unsure
# Is it a gpcr based on title?
sub title_type {
	my $t = shift;

	# search negative first
	foreach $p (@title_no) {
		if ($t =~ /$p/i) { return 0; }
	}

	# now the positive
	foreach $p (@title_yes) {
		if ($t =~ /$p/i) { return 1; }
	}

	if ($t =~ /Hypothetical protein/i) {
	# just don't want to do anything
		return 2;
	}

	return 3;
}

sub keyword_type {
	my $kwd_r = shift;
	foreach $k (@$kwd_r) {
		foreach $p (@keyword_no) {
			if ($k =~ /$p/i) { return 0; }
		}
		if ($k =~ /G.protein.coupled receptor/i) { return 1; }
	}
	return 3;
}

sub toarr {
	my $ln = shift;
	chomp $ln;          # remove \n
	$ln =~ s/\";$|\"$//;
	$ln = substr($ln,1); # remove first \"
	return split(/$delimiter/, $ln);
}


=for fromPGDB
	This method is very slow.
=cut
sub fromPGDB {
	## first we need to get the input candidates
	$pgdb = get_orparadb;
	$result = $pgdb->exec("select distinct target from stdblp_gpcr where valid isnull");
	if ($result->resultStatus != PGRES_TUPLES_OK) {
		die "query failed", $pgdb->errorMessage, "\n";
	}
	$swdb = get_acedb("swprt");

	$i=0;
	while ($i < $result->ntuples) {
		$target=$result->getvalue($i,0);
		$prt = $swdb->fetch(protein=>$target) or die "$target not found in SWISSPRT\n";

		print "protein $target\n";
		print $target, " ", $prt->get("title")->at, "\n";
		$i++;
	}

}
