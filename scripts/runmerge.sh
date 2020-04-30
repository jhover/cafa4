#!/bin/bash
#
# create /allspecies.noknow.2010.allmethods.allaspects.eval.csv
#
INDIR=~/play/cafa4/ismb2
OUTDIR=~/play/cafa4/ismb2
PROG=~/data/git/cafa4/scripts/mergeevals.py
KNOWLEDGE="noknow limited"

echo "running merge..."

for KNOW in $KNOWLEDGE; do 
	echo "Handling knowledge $KNOW..."
	OUTFILE=allspecies.$KNOW.2010.allmethods.allaspects.eval.csv
	echo "$PROG $INDIR/$KNOW/sp_species.*{bp,mf,cc}.eval.csv > $OUTDIR/$OUTFILE "
	$PROG $INDIR/$KNOW/sp_species.*{bp,mf,cc}.eval.csv > $OUTDIR/$OUTFILE 

done
