#!/bin/bash
if [ "$#" -ne 2 ]; then
        echo "Please specify 1) input folder containing indexeded bam files 2) output folder for bigwig files"
	exit 1
fi

inDir=$1
outDir=$2

mkdir -p ${outDir}

echo "Generate bigWig files using bamCoverage of deepTools package."
for FILE in ${inDir}/*"_dedup.bam";
do
    filename=${FILE##*/}
    filename=${filename/"_dedup.bam"/}
    output=${outDir}/$filename"_dedup.bigWig"
    
    
    if [ -f "$output" ];
    then
        echo "$output already exists"
        continue
    fi
    
    if [[ "$FILE" = *"_dedup.bigWig" ]];
    then
        echo "$FILE already exists"
        continue
    fi
    
    bamCoverage --normalizeUsing RPKM --ignoreDuplicates -e 0 -bs 10 -b $FILE -o $output"

done
