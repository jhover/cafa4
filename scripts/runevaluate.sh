#!/bin/bash
#
#  Run evaluations on all predictions in testdir
#  by aspect and all. 
#
#INDIR=~/play/ismb/noknow
OUTDIR=~/play/ismb/evaluate
PROG=~/data/git/cafa4/fastcafa/fastcafa.py
CONF=~/data/git/cafa4/etc/fastcafa.conf

KNOWLEDGE="noknow"
METHODS="expression phmmer orthoexpression prior"
ASPECTS="bp mf cc all"

#VERSION="2010"
VFLAG=" -V $VERSION "
#DEBUG=" -d "
DEBUG=" "

mkdir -p $OUTDIR
echo "Running on all targets..."

for KNOW in $KNOWLEDGE; do
	echo "Handling knowledge $KNOW..."
	OUTBASE=$OUTDIR/$KNOW
	echo "Making outdir $OUTBASE..."
	mkdir -p $OUTBASE
	for METHOD in $METHODS; do 
		echo "Handling method $METHOD..."
		
		for ASPECT in $ASPECTS; do
			echo "Handling aspect $ASPECT..." 
		
			for FILE in `ls $INDIR/*.$KNOW.*.$METHOD.$ASPECT.csv`; do 
				echo "Handling $FILE..."
				echo "###############################################"
				FILENAME=`basename $FILE`
				EXTENSION="${FILENAME##*.}"
				FILEBASE="${FILENAME%.*}"
				OUTFILE=$OUTBASE/$FILEBASE.eval.csv
				echo "FILE=$FILE"
				echo "EXTENSION=$EXTENSION"
				echo "FILEBASE=$FILEBASE"
				echo "OUTBASE=$OUTBASE"
				echo "OUT=$OUTFILE"
				if [ ! -f $OUTFILE ]; then
					echo "time $PROG -C $DEBUG -c $CONF evaluate -g $ASPECT -i $FILE -o $OUTFILE "
					echo ""
					time $PROG -C $DEBUG -c $CONF evaluate -g $ASPECT -i $FILE -o $OUTFILE
				else
					echo "Output exists. Skipping..."				
				fi
			done
		done
	done
done	
