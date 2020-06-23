#!/bin/bash
#
#
# SPECLIST="HUMAN MOUSE RAT BOVIN PIG DANRE DROME CAEEL ARATH MAIZE YEAST SCHPO ECOLI PSEAI  DICDI MYCGE BACCR SALTY"
#SPECLIST="ARATH CAEEL DANRE DROME HUMAN RAT YEAST"
SPECLIST="ARATH CAEEL DANRE DROME HUMAN MOUSE RAT YEAST"
PROG=~/git/cafa4/fastcafa/fastcafa.py
CONF=~/git/cafa4/etc/fastcafa.conf
NUMSEQ=300
OLD=2010
NEW=2019
OUTDIR=~/play/ismb/testtargets/$OLD-$NEW

mkdir -p $OUTDIR

for SPECIES in $SPECLIST; do 
	TFAOUT=$OUTDIR/sp_species.$SPECIES.test.noknow.$NUMSEQ.tfa
	echo "time $PROG -c $CONF testset -n $NUMSEQ -s $SPECIES -o $TFAOUT -O $OLD -N $NEW  "  
	time $PROG -c $CONF testset -n $NUMSEQ -s $SPECIES -o $TFAOUT -O $OLD -N $NEW 
done

for SPECIES in $SPECLIST; do 
	TFAOUT=$OUTDIR/sp_species.$SPECIES.test.limited.$NUMSEQ.tfa
	echo "time $PROG -c $CONF testset -n $NUMSEQ -s $SPECIES -L -o $TFAOUT -O $OLD -N $NEW  "  
	time $PROG -c $CONF testset -n $NUMSEQ -s $SPECIES -L -o $TFAOUT -O $OLD -N $NEW 
done

