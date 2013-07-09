#!/usr/bin/sh

/package/UCSCTools/bigWigToBedGraph $1 _temp.bed
cat _temp.bed | perl -lane '$id+=1; if($F[3]<1) { print "$F[0]\t$F[1]\t$F[2]\tunmap_$id\t0"}' |  /package/BEDTools/bin/mergeBed > $2
rm _temp.bed
