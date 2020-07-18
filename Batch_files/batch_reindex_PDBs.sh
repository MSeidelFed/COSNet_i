#!/bin/bash
#Batch_reindex_PDBs.sh: 

PYCODEPATH=$1
FILENAME=$2
PDBPATH=$3
OUTPATH=$4
STREIXNR=$5
COUNT=0
IFS=" "

[ $# -eq 0 ] && { echo "-----------------------------------------------------------------------------------------"; echo "Usage: bash $0 [pathtopythoncode] [pathto/infile] [pathtopdbs] [outpath] [startingreindexnumber]"; echo "-----------------------------------------------------------------------------------------"; exit 1; }

while read line
do
	COUNT=`expr $COUNT + 1`
	echo "Running on protein # "$COUNT 
	read -ra ADDR <<< $line
	echo "${ADDR[0]}"
	echo "${ADDR[1]}"
	python3 -u "$PYCODEPATH"reindex_pdb.py $STREIXNR $PDBPATH${ADDR[0]} $OUTPATH${ADDR[1]}  

done < $FILENAME
