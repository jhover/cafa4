#!/bin/bash
#
# The CAFA target species....
# Should normally take about 20 minutes to run. 
#
SPECLIST="ARATH BACCR BOVIN CAEEL DANRE DICDI DROME ECOLI HUMAN MAIZE MOUSE MYCGE PIG PSEAI SCHPO RAT SALTY YEAST "
VERSIONS="current previous"
PROG=~/data/git/cafa4/fastcafa/fastcafa.py
CONF=~/data/git/cafa4/etc/fastcafa.conf
OUTDIR=~/play/cafa4

for VERSION in $VERSIONS; do 
	echo "Building prior(s) for uniprot version $VERSION..."
	for SPECIES in $SPECLIST; do 
		echo "Building prior for $SPECIES..." 
		echo "time $PROG -c $CONF build_prior -s $SPECIES -V $VERSION -o $OUTDIR/uniprot.prior.$VERSION.$SPECIES.csv" 
		time $PROG -c $CONF build_prior -s $SPECIES -V $VERSION -o $OUTDIR/uniprot.prior.$VERSION.$SPECIES.csv    
	done
	
	echo "Building global prior..."
	echo "time $PROG -c $CONF build_prior -V $VERSION -o $OUTDIR/uniprot.prior.$VERSION.GLOBAL.csv"
	time $PROG -c $CONF build_prior -V $VERSION -o $OUTDIR/uniprot.prior.$VERSION.GLOBAL.csv

done
