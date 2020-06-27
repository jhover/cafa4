#!/bin/bash
#
# Creates symlinks in species subdirs in order to 
# run all species in parallel
# Exports INDIR to env for runevaluate.sh to pick it up. 
#

WORKDIR=~/play/ismb
INDIR=$WORKDIR/predict

SPECLIST="ARATH CAEEL DANRE DROME HUMAN MOUSE RAT YEAST"

cd $WORKDIR
for SPECIES in $SPECLIST; do
    mkdir -p $SPECIES
    cd $SPECIES
    echo ln -s $INDIR/*.$SPECIES.* ./
    ln -s $INDIR/*.$SPECIES.* ./
    cd ..

    #echo mv -v $INDIR/*.$SPECIES.* $WORKDIR/$SPECIES/
    # mv -v $INDIR/*.$SPECIES.* $WORKDIR/$SPECIES/
done

cd $WORKDIR

for SPECIES in $SPECLIST; do
	export INDIR=~/play/ismb/$SPECIES
	echo $INDIR
	echo "time ~/git/cafa4/scripts/runevaluate.sh > runevaluate.${SPECIES}.log 2>&1 &"
	time ~/git/cafa4/scripts/runevaluate.sh	> runevaluate.${SPECIES}.log 2>&1 &
done

#for SPECIES in $SPECLIST; do
#    echo mv -v $WORKDIR/$SPECIES/*.${SPECIES}.* $WORKDIR/predict
#    mv -v $WORKDIR/$SPECIES/*.${SPECIES}.* $WORKDIR/predict
# done
