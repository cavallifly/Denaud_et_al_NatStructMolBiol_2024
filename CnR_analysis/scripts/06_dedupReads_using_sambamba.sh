#!/bin/bash
if [ "$#" -ne 2 ];
then
    echo "Please specify 1) input folder containing bam files and 2) output folder for text files"
    exit 1
fi

inDir=$1
outDir=$2

mkdir -p ${outDir}

echo "Removing duplicated reads using sambamba"
for FILE in ${inDir}/*"_s.bam"; do
	filename=${FILE##*/}
	filename=${filename/".bam"/}
	output=${outDir}/${filename}"_dedup.bam"

	if [ -f "$output" ]; then
                echo "$output already exists"
                continue
        fi
     
	sambamba markdup -r --hash-table-size 500000 --overflow-list-size 500000 $FILE $output"

done
