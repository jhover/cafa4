#!/bin/bash
#
# Run combine combinations on basenames in provided file. 
#
#  E.g. ls | awk -F. '{ print $1"."$2"."$3"."$4 }' | sort | uniq > infile
# 

BASEFILE=$1

PROG=~/git/cafa4/fastcafa/fastcafa.py
CONF=~/git/cafa4/etc/fastcafa.conf
 
DEBUG="-d "
#DEBUG=""

BASES=`cat $1`

for FILEBASE in $BASES; do 

	echo "Handling $FILEBASE..."
		
	echo  "###############################################"
	echo "Running combines..."
	echo "Doing phmmer_prior..."
	INONE=$FILEBASE.phmmer.csv
	INTWO=$FILEBASE.prior.csv
	OUTFILE=$FILEBASE.phmmer_prior.csv
	if [ ! -f $OUTFILE ]; then 
		echo "time $PROG -C $DEBUG -c $CONF combine -f $INONE -s $INTWO -o $OUTFILE  "
		echo ""
		time $PROG -C $DEBUG -c $CONF combine -f $INONE -s $INTWO -o $OUTFILE 
	fi
	
	echo "Doing phmmer_expression..."
	INONE=$FILEBASE.phmmer.csv
	INTWO=$FILEBASE.expression.csv
	OUTFILE=$FILEBASE.phmmer_expression.csv
	if [ ! -f $OUTFILE ]; then 
		echo "time $PROG -C $DEBUG -c $CONF combine -f $INONE -s $INTWO -o $OUTFILE  "
		echo ""
		time $PROG -C $DEBUG -c $CONF combine -f $INONE -s $INTWO -o $OUTFILE  
	fi
	
	echo "Doing phmmer_expression_prior..."
	INONE=$FILEBASE.phmmer_expression.csv
	INTWO=$FILEBASE.prior.csv
	OUTFILE=$FILEBASE.phmmer_expression_prior.csv
	if [ ! -f $OUTFILE ]; then
		echo "time $PROG -C $DEBUG -c $CONF combine -f $INONE -s $INTWO -o $OUTFILE "
		echo ""
		time $PROG -C $DEBUG -c $CONF combine -f $INONE -s $INTWO -o $OUTFILE  
	fi
	
	echo "Doing phmmer_prior_expression..."
	INONE=$FILEBASE.phmmer_prior.csv
	INTWO=$FILEBASE.expression.csv
	OUTFILE=$FILEBASE.phmmer_prior_expression.csv
	if [ ! -f $OUTFILE ]; then
		echo "time $PROG -C $DEBUG -c $CONF combine -f $INONE -s $INTWO -o $OUTFILE "
		echo ""
		time $PROG -C $DEBUG -c $CONF combine -f $INONE -s $INTWO -o $OUTFILE  
	fi
	
	echo "Doing expression_phmmer..."
	INONE=$FILEBASE.expression.csv
	INTWO=$FILEBASE.phmmer.csv
	OUTFILE=$FILEBASE.expression_phmmer.csv
	if [ ! -f $OUTFILE ]; then
		echo "time $PROG -C $DEBUG -c $CONF combine -f $INONE -s $INTWO -o $OUTFILE  "
		echo ""
		time $PROG -C $DEBUG -c $CONF combine -f $INONE -s $INTWO -o $OUTFILE  
	fi

done

