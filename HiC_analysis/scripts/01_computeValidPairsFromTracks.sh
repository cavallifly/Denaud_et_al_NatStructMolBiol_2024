inFile=01_computeValidPairsFromTracks_genomeWide.out
Rscript ./scripts/01_computeValidPairsFromTracks.R > ${inFile}

outFile=count_validPairs_genomeWide.out
grep hic_ ${inFile} | grep -v "chr\|dac" | awk '{if(NF==3) print $0}' | awk '{printf("%s\t%s\t%s\t%s\n",$2,$3,$4,$5)}' | sed "s/\"//g" > ${outFile}
echo
outFile=count_validPairs_chr2L.out
grep hic_ ${inFile} | grep "chr2L chr2L" | awk '{printf("%s\t%s\t%s\t%s\n",$2,$3,$4,$5)}' | sed "s/\"//g" > ${outFile}
echo
outFile=count_validPairs_dac.out
grep hic_ ${inFile} | grep dac           | awk '{printf("%s\t%s\t%s\t%s\n",$2,$3,$4,$5)}' | sed "s/\"//g" > ${outFile}
