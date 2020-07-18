#!/bin/bash
#batch_calc_dist.sh: 

PYCODEPATH=$1
FILENAME=$2
PDBPATH=$3
OUTPATH=$4
COUNT=0
IFS=" "

[ $# -eq 0 ] && { echo "-----------------------------------------------------------------------------------------"; echo "Usage: bash $0 [pathtopythoncode] [pathto/infile] [pathtopdbs] [outpath]"; echo "-----------------------------------------------------------------------------------------"; exit 1; }

while read line
do
	COUNT=`expr $COUNT + 1`
	echo "Running on protein # "$COUNT 
	read -ra ADDR <<< $line
	echo "${ADDR[0]}"
	echo "${ADDR[1]}"
	python3 -u "$PYCODEPATH"calculate_distance.py $PDBPATH${ADDR[0]} $PDBPATH${ADDR[1]} $OUTPATH  

done < $FILENAME
