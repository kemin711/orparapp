#!/usr/bin/perl -w

#convert the stupid GFF format back into ace

$_ = <>;
while (!/^##sequence-region/) {$_=<>;}
@arr=split / /;
$seq = $arr[1];
$begin = $arr[2];
$end = $arr[3];

$_ = <>;  #summary line not needed!!
$_ = <>;  #first real data
@arr = split /\t/;
#################################collect data
$seqName = $arr[0];
$evidence = $arr[1];
$featureType = $arr[2];
$featBegin = $arr[3];
$featEnd = $arr[4];
$score = $arr[5];
$strand = $arr[6];
$frame = $arr[7];
$attribute = $arr[8];
#############################

print "position of this sequence in Chr22 $begin   $end\n";
if ($seqName ne $seq) { 
	print "something may be wrong with sequence name\n";;
}


print "Sequence $seq\n";
print "Evidence $evidence\n";
if ($strand eq "+") {
	print "Subsequence $featureType $featBegin $featEnd\n";
}
else {
	print "Subsequence $featureType $featEnd $featBegin\n";
}
print "Title $attribute\n";





