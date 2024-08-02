#!/bin/bash

bam1=$1
bam2=$2



output=mergedBam.bam
	
samtools merge -o $output $bam1 $bam2 $bam3 $bam4
samtools index ${output}
