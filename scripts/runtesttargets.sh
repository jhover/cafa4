#!/bin/bash
#
#  Make predictions(prior,  phmmer and expression) on input files using previous 
#  Output .csv predictions. 
#
TESTDIR=~/data/cafa4/testtargets
OUTDIR=~/play/cafa4/evaluate
PROG=~/data/git/cafa4/fastcafa/fastcafa.py
CONF=~/data/git/cafa4/etc/fastcafa.conf
METHODS="prior phmmer expression"
#DEBUG=" -d "
DEBUG=""
VERSION="previous"
VFLAG=" -V $VERSION "

echo "Running on all targets..."
for TFA in `ls $TESTDIR/*.tfa`; do
	echo "Handling $TFA..."
	echo "###############################################"
	FILENAME=`basename $TFA`
	EXTENSION="${FILENAME##*.}"
	FILEBASE="${FILENAME%.*}"
	#echo "$FILENAME $FILEBASE $EXTENSION"
	for METHOD in $METHODS; do
		# echo "method is $METHOD"
		PREDOUT=$OUTDIR/$FILEBASE.$VERSION.$METHOD.csv
		
		echo "running $METHOD ..."
		echo "time $PROG -C $DEBUG -c $CONF $METHOD $VFLAG -i $TFA -o $PREDOUT  "
		echo ""
		time $PROG -C $DEBUG -c $CONF $METHOD $VFLAG -i $TFA -o $PREDOUT 

	done
	echo  "###############################################"
done
