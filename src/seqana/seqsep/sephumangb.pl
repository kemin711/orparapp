### perl ###
##!/bin/sh
# file: sephumangb
# you have to call this from the directory containing all the gbpri files

files=`ls gbpri*.seq`
for f in $files
	do echo separating human sequence from $f ....
		seqsep $f -org human -s 5000
		echo
	done

