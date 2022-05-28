#!/usr/bin/env bash

#Checking for dependencies

if ! command -v minimap2 > /dev/null; then
        printf "\n Minimap2 not found in path. Exiting"
        exit
fi

if ! command -v bwa > /dev/null; then
	printf "\n BWA not found in path. Exiting"
	exit
fi

if ! command -v racon > /dev/null; then
	printf "\n Racon not found in path. Exiting"
	exit
fi

if ! command -v medaka_consensus > /dev/null; then
	printf "\n Medaka not found in path. Exiting"
	exit
fi

if ! command -v polypolish > /dev/null; then
	printf "\n Polypolish not found in path. Exiting"
	exit
fi


#Checking to see if there are fasta files in the directory

count=`ls -1 *.fasta 2>/dev/null | wc -l`
if [ $count != 0 ]
then
	echo "Fasta files found - processing"
else
	echo "No fasta files found in directory. Exiting"
	exit
fi

threads=$1
long_read=$(echo $2 | tr -d '\r')
short_read_1=$(echo $3 | tr -d '\r')
short_read_2=$(echo $4 | tr -d '\r')

for F in *.fasta; do
	N=$(basename $F .fasta);
	minimap2 -t "$threads" -ax map-ont $F "$long_read" > $N.sam;
	racon -m 8 -x -6 -g -8 -w 500 -t 10 "$long_read" $N.sam $F > racon_$N.fna;
	medaka_consensus -i "$long_read" -d racon_$N.fna -o medaka_$N -t "$threads" -m r941_min_sup_g507;
	rm $N.sam racon_$N.*;
	if [[ $# == 4 ]]; then
		echo "#####################################################################";
		echo "                                                                     ";
		echo "Short reads detected as arguments 2 and 3. Proceeding with polypolish";
		echo "                                                                     ";
		echo "#####################################################################";
		mkdir polypolish_$N;
		bwa index medaka_$N/consensus.fasta;
		bwa mem -t "$threads" -a medaka_$N/consensus.fasta "$short_read_1" > polypolish_$N/alignment1.sam;
		bwa mem -t "$threads" -a medaka_$N/consensus.fasta "$short_read_2" > polypolish_$N/alignment2.sam;
		polypolish medaka_$N/consensus.fasta polypolish_$N/alignment1.sam polypolish_$N/alignment2.sam > polypolish_$N/polypolish_$F;
		rm polypolish_$N/alignment1.sam polypolish_$N/alignment2.sam;
		rm -rf medaka_$N;
	else
		echo "Short reads not provided for $F";
	fi;
done

	echo "#############################"
	echo "                             "
	echo "$F polishing finished"
	echo "                             "
	echo "#############################"

