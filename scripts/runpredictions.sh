#!/bin/bash
#
#  Make predictions on input files using old version of uniprot.
#  Use general + all aspects individually.   
#  Output .csv predictions. 
#
# SPECIES="BOVIN DICDI ECOLI"
#SPECIES="DICDI"


TESTDIR=~/play/ismb/testtargets/2010-2019
OUTDIR=~/play/ismb/
PROG=~/data/git/cafa4/fastcafa/fastcafa.py
CONF=~/data/git/cafa4/etc/fastcafa.conf
#METHODS="prior phmmer expression orthoexpression"
METHODS="prior phmmer expression orthoexpression"
ASPECTS="bp cc mf all"

#DEBUG=" -d "
DEBUG=" "
VERSION="2010"
VFLAG=" -V $VERSION "

mkdir -p $OUTDIR
echo "Running on all targets..."
for TFA in `ls $TESTDIR/*.${SPECIES}.*.tfa`; do
	echo "Handling $TFA..."
	echo "###############################################"
	FILENAME=`basename $TFA`
	EXTENSION="${FILENAME##*.}"
	FILEBASE="${FILENAME%.*}"
	for METHOD in $METHODS; do
		echo "Handling method $METHOD for all aspects..."

		for ASPECT in $ASPECTS; do
			echo "Handling aspect $ASPECT..."
			PREDOUT=$OUTDIR/$FILEBASE.$VERSION.$METHOD.$ASPECT.csv		
			if [ -f "$PREDOUT" ]; then
				echo "Output exists. Skipping..."
			else
				echo "Handling method $METHOD with aspect $ASPECT..."
				echo "time $PROG -C $DEBUG -c $CONF $METHOD -g $ASPECT $VFLAG -i $TFA -o $PREDOUT  "
				echo ""
				time $PROG -C $DEBUG -c $CONF $METHOD -g $ASPECT $VFLAG -i $TFA -o $PREDOUT 
			fi
		done
	done
	echo  "###############################################"
done
