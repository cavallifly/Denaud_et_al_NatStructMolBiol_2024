#!/bin/bash
if [ "$#" -ne 2 ];
then
    echo "Please specify 1) input folder containing fq.gz files to do fastQC 2)output folder"
    exit 1
fi

inDir=$1"/"
outDir=$2"/"

mkdir -p ${outDir}

echo "Doing fastQC analysis on the files in folder $1."
echo "The output will be stored in the folder $2"
for FILE in ${inDir}/*"fq.gz";
do
    echo ${FILE##*/}

    bash fastqc $FILE --outdir=${outDir}

done
