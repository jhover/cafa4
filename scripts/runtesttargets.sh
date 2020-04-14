#!/bin/bash
#
#  Make predictions(prior,  phmmer and expression) on input files using previous 
#  Output .csv predictions. 
#
TESTDIR=~/data/cafa4/testtargets/2017-2019/
OUTDIR=~/play/cafa4/20200408
PROG=~/data/git/cafa4/fastcafa/fastcafa.py
CONF=~/data/git/cafa4/etc/fastcafa.conf
METHODS="expression"
DEBUG=" -d "
#DEBUG=""
VERSION="2017"
VFLAG=" -V $VERSION "

mkdir -p $OUTDIR


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
