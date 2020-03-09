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
	PREDOUT=~/play/cafa4/sp_species.$SPECIES.test.$NUMSEQ.expr.previous.csv
	EVALOUT=~/play/cafa4/sp_species.$SPECIES.test.$NUMSEQ.evaluate.csv

	echo "time $PROG -c $CONF testset -n $NUMSEQ -s $SPECIES -o $TFAOUT"  
	time $PROG -c $CONF testset -n $NUMSEQ -s $SPECIES -o $TFAOUT

	echo "time $PROG -C -c $CONF expression -i $TFAOUT -o $PREDOUT -V previous"
	time $PROG -C -c $CONF expression -i $TFAOUT -o $PREDOUT -V previous

	echo "time $PROG -C -c $CONF evaluate -p $PREDOUT -o $EVALOUT"
	time $PROG -C -c $CONF evaluate -p $PREDOUT -o $EVALOUT
done


