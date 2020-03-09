#!/bin/bash
#
#  Run eval (phmmer and expression) on input files. 
#  Output .csv predictions. 
#
TESTDIR=~/data/cafa4/testtargets
OUTDIR=~/play/cafa4/evaluate
PROG=~/git/cafa4/fastcafa/fastcafa.py
CONF=~/git/cafa4/etc/fastcafa.conf
METHODS="prior phmmer expression"
DEBUG="-d "


echo "Running tests...."

for TFA in `ls $TESTDIR/*.tfa`; do
	echo "Handling $TFA..."
	echo "###############################################"
	FILENAME=`basename $TFA`
	EXTENSION="${FILENAME##*.}"
	FILEBASE="${FILENAME%.*}"
	#echo "$FILENAME $FILEBASE $EXTENSION"
	for METHOD in $METHODS; do
		# echo "method is $METHOD"
		PREDOUT=$OUTDIR/$FILEBASE.$METHOD.csv
		EVALOUT=$OUTDIR/$FILEBASE.$METHOD.eval.csv
		
		echo "running $METHOD ..."
		echo "time $PROG -C $DEBUG -c $CONF $METHOD -i $TFA -o $PREDOUT -V previous"
		echo ""
		time $PROG -C $DEBUG -c $CONF $METHOD -i $TFA -o $PREDOUT -V previous
				
		echo "Running evaluate..."
		if [ -f $PREDOUT ]; then
			echo "time $PROG -C $DEBUG -c $CONF evaluate -p $PREDOUT -o $EVALOUT"
			echo ""
			time $PROG -C $DEBUG -c $CONF evaluate -p $PREDOUT -o $EVALOUT
		else
			echo "no infile $PREDOUT skipping..."
		fi
	done
	echo  "###############################################"
done

 

 