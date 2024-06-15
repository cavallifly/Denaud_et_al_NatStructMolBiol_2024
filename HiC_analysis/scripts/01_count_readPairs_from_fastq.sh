outFile=${PWD}/count_readPairs_from_fastq.out
for dir in $(ls -1 | grep hic_);
do    
    cd $dir
    echo $dir
    for file in $(ls -1 *fq 2> /dev/null);
    do
	echo $file
	nLines=$(grep $file ${outFile} | wc -l)
	echo $nLines
	if [[ ${nLines} -ne 0 ]];	   
	then
	    continue
	fi
	echo $file
	echo "$dir $file Nreads $(wc -l $file | awk '{print $1/4}')" >> ${outFile}	
    done
    for file in $(ls -1 *gz 2> /dev/null);
    do
	echo $file
	nLines=$(grep $file ${outFile} | wc -l)
	echo $nLines
	if [[ ${nLines} -ne 0 ]];	   
	then
	    continue
	fi
	echo $file
	echo "$dir $file Nreads $(zcat $file | wc -l | awk '{print $1/4}')" >> ${outFile}
    done    
    cd ..
done
