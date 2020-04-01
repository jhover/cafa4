#!/bin/bash
#
#  Run eval (phmmer and expression) on input files. 
#  Output .csv predictions. 
#

TFA=$1
OUTDIR=$3


# RAT
#TFA="$TESTDIR/sp_species.10116.bmtest.limited.20.tfa"
#TFA="$TESTDIR/sp_species.10116.bmtest.noknow.2.tfa"

# MOUSE
#TFA="$TESTDIR/sp_species.10090.bmtest.noknow.20.tfa"
#TFA="$TESTDIR/sp_species.10090.bmtest.limited.20.tfa"

# OUTDIR=~/play/cafa4

PROG=~/git/cafa4/fastcafa/fastcafa.py
CONF=~/git/cafa4/etc/fastcafa.conf
METHODS="prior phmmer expression"
PREDTYPES="prior phmmer expression phmmer_prior phmmer_expression phmmer_expression_prior phmmer_prior_expression expression_phmmer"
# PREDTYPES="prior phmmer expression phmmer_prior phmmer_prior_expression"  
#DEBUG="-d "
DEBUG=""

echo "Running full pipeline on on file..."
echo "Handling $TFA..."
echo "###############################################"
FILENAME=`basename $TFA`
EXTENSION="${FILENAME##*.}"
FILEBASE="${FILENAME%.*}"
#echo "$FILENAME $FILEBASE $EXTENSION"

if [ ! -f $TFA ]; then
	echo "file not found"
	exit 1
fi


echo "################################################"
echo "Running evaluates..."

for PREDTYPE in $PREDTYPES ; do 
	echo "Evaluating prediction $PREDTYPE..."
	INPRED=$OUTDIR/$FILEBASE.$PREDTYPE.csv
	EVALOUT=$OUTDIR/$FILEBASE.$PREDTYPE.eval.csv
	if [ ! -f $EVALOUT ]; then
		echo "time $PROG -C $DEBUG -c $CONF evaluate -i $INPRED -o $EVALOUT "
		echo ""
		time $PROG -C $DEBUG -c $CONF evaluate -i $INPRED -o $EVALOUT
	fi
done

echo "################################################"
# echo "Outputting CAFA file..."
