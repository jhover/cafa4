#!/bin/bash
#
#
SPECLIST="HUMAN MOUSE RAT BOVIN PIG DANRE DROME CAEEL ARATH MAIZE YEAST SCHPO ECOLI PSEAI  DICDI MYCGE BACCR SALTY"
#SPECLIST="ARATH CAEEL DANRE DROME HUMAN RAT YEAST"
PROG=~/git/cafa4/fastcafa/fastcafa.py
CONF=~/git/cafa4/etc/fastcafa.conf
NUMSEQ=2


for SPECIES in $SPECLIST; do 
	TFAOUT=~/play/cafa4/sp_species.$SPECIES.test.$NUMSEQ.tfa

	echo "time $PROG -c $CONF testset -n $NUMSEQ -s $SPECIES -o $TFAOUT"  
	time $PROG -c $CONF testset -n $NUMSEQ -s $SPECIES -o $TFAOUT


done


