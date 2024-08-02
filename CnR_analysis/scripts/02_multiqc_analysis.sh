#!/bin/bash
if [ "$#" -ne 2 ]; then
        echo "Please specify 1) input folder containing fastqc.zip files to do multiQC 2)output folder"
        exit 1
fi

inDir=$1
outDir=$2

echo "Doing multiQC analysis on the files in $1."
echo "The output will be stored in $2."

mkdir -p ${outDir}

multiqc ${inDir}/*fastqc.zip -o ${outDir} --pdf

