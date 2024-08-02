#!/bin/bash
if [ "$#" -ne 2 ];
then
    echo "Please specify 1) input folder containing sam files and 2) output folder for mapq-sorted .bam files"
    exit 1
fi

inDir=$1
outDir=$2

mkdir -p ${outDir}

echo "Filtering out reads with MAPQ < 30."
for FILE in $(ls -1 ${inDir}/*.sam);
do
    filename=${FILE##*/}
    filename=${filename/".sam"/}
    output=${outDir}/$filename".bam"
    
    if [ -f "$output" ];
    then
        echo "$output already exists"
        continue
    fi

    echo "$cmd"
    bash samtools view -b -q 30 $FILE -o $output
done
