#!/bin/bash
if [ "$#" -ne 3 ];
then
    echo "Please specify"
    echo "1) input folder containing fq.gz files"
    echo "2) SAM output folder"
    echo "3) Bowtie2Index folder for dm6 assembly"
    exit 1
fi

inDir=$1
outDir=$2
bowtie2index=$3


echo "Doing bowtie2 alignment on dm6 for the samples in ${inDir}."
for FILE in ${inDir}/*_1.fq.gz;
do

    fq=${FILE##*/}
    name=${fq/".fq.gz"/}
    sam=${outDir}/$name".sam"

    if [ -f "$sam" ];
    then
        echo "$sam already exists"
        continue
    fi

    report=${outDir}/${name/_1/}.bowtiereport.txt
    
    if [[ "$FILE" = *"_1.fq.gz" ]];
    then
	cmd="bowtie2 -p 40 -x ${bowtie2index} --local --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700 -1 "$FILE" -2 "${FILE/_1.fq.gz/_2.fq.gz}" -S "$sam" 2>"$report""
    else
	cmd="bowtie2 -p 12 -x ${bowtie2index} --local --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700 -U "$FILE" -S "$sam" >> "$report" 2>&1"
    fi

    echo "Running command $cmd"
    bash $cmd

done
