#!/bin/bash
PROG=~/git/cafa4/fastcafa/fastcafa.py
CONF=~/git/cafa4/etc/fastcafa.conf
DEBUG="-d "

FILES=$@

for FILE in $FILES; do 
	echo "Handling $FILE..."
	echo "###############################################"
	FILENAME=`basename $FILE`
	EXTENSION="${FILENAME##*.}"
	FILEBASE="${FILENAME%.*}"
	OUTFILE=$FILEBASE.eval.csv
	echo "FILE=$FILE"
	echo "EXTENSION=$EXTENSION"
	echo "FILEBASE=$FILEBASE"
	echo "OUT=$OUTFILE"
	
	if [ ! -f $OUTFILE ]; then
		echo "time $PROG -C $DEBUG -c $CONF evaluate -i $FILE -o $OUTFILE "
		echo ""
		time $PROG -C $DEBUG -c $CONF evaluate -i $FILE -o $OUTFILE
	fi


done
