### perl ###
##!/bin/sh

# start with gbvrt.seq file

echo separating the sau division from all vertebrates ...
echo Sauropsids
seqdiv -s 10000 -f gbvrt.seq -d sau Archosauria Lepidosauria Testudines
mv gbvrt.sau sau
mv gbvrt.rest gbvrt

echo separating the amp division from all vertebrates ...
echo Amphibia separation
seqdiv -s 8000 -f gbvrt -d amp Amphibia
mv gbvrt.amp amp
#rm gbvrt
mv gbvrt.rest gbvrt
seqsep -s 10000 -org Amphibia gbvrt  # two entries with Organism=Amphibia
mv gbvrt.amphibia amp/amp_poor.seq
#rm gbvrt
mv gbvrt.other gbvrt
echo

echo removing fugu from all vertebrates ...
seqsep -s 8000 -org fugu gbvrt
rm gbvrt.fugu
mv gbvrt.other gbvrt
echo


# needs to remove all organisms higher than amp to be separated
# to carry out this step
echo separating the bony division from all vertebrates ...
echo bony fishes separation
seqdiv -s 8000 -f gbvrt -d bony Actinopterygii Coelacanthiformes Dipnoi # Euteleostomi
mv gbvrt.bony bony
#rm gbvrt
mv gbvrt.rest gbvrt
echo

echo separating the cart division from all vertebrates ...
echo cartilaginous fishes sep
seqdiv -f gbvrt -d cart Chondrichthyes
mv gbvrt.cart cart
echo
