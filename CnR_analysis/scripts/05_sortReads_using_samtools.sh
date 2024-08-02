#!/bin/bash
if [ "$#" -ne 2 ];
then
    echo "Please specify 1) input folder containing bam files and 2) output folder for sorted _s.bam files"
    exit 1
fi

inDir=$1
outDir=$2

mkdir -p ${outDir}

echo "sort bam files"
for FILE in ${inDir}/*".bam";
do
    filename=${FILE##*/}
    filename=${filename/".bam"/}
    output=${outDir}/$filename"_s.bam"
    
    if [ -f "$output" ];
    then
        echo "$output already exists"
        continue
    fi
    
    if [[ "$FILE" = *"_s.bam" ]];
    then
        echo "$FILE already exists"
        continue
    fi
    
    cmd="samtools sort $FILE -o $output"
    echo "Running command $cmd"
    bash $cmd

done
