#!/bin/bash
#
#  Make predictions on input files using old version of uniprot.
#  Use general + all aspects individually.   
#  Output .csv predictions. 
#
TESTDIR=~/data/cafa4/testtargets/2010-2019
OUTDIR=~/play/cafa4/20200414
PROG=~/data/git/cafa4/fastcafa/fastcafa.py
CONF=~/data/git/cafa4/etc/fastcafa.conf
METHODS="prior"
ASPECTS="bp cc mf"
#DEBUG=" -d "
DEBUG=" "
VERSION="2010"
VFLAG=" -V $VERSION "

mkdir -p $OUTDIR
echo "Running on all targets..."
for TFA in `ls $TESTDIR/*.tfa`; do
	echo "Handling $TFA..."
	echo "###############################################"
	FILENAME=`basename $TFA`
	EXTENSION="${FILENAME##*.}"
	FILEBASE="${FILENAME%.*}"
	for METHOD in $METHODS; do
		echo "Handling method $METHOD for all aspects..."
		PREDOUT=$OUTDIR/$FILEBASE.$VERSION.$METHOD.all.csv		
		echo "running $METHOD ..."
		echo "time $PROG -C $DEBUG -c $CONF $METHOD $VFLAG -i $TFA -o $PREDOUT  "
		echo ""
		time $PROG -C $DEBUG -c $CONF $METHOD $VFLAG -i $TFA -o $PREDOUT 

		for ASPECT in $ASPECTS; do
			PREDOUT=$OUTDIR/$FILEBASE.$VERSION.$METHOD.$ASPECT.csv		
			echo "Handling method $METHOD with aspect $ASPECT..."
			echo "time $PROG -C $DEBUG -c $CONF $METHOD -g $ASPECT $VFLAG -i $TFA -o $PREDOUT  "
			echo ""
			time $PROG -C $DEBUG -c $CONF $METHOD -g $ASPECT $VFLAG -i $TFA -o $PREDOUT 
		done
	done
	echo  "###############################################"
done
