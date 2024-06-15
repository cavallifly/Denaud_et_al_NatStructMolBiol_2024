size=800

for inFile in $(ls -1 *toModify*.png);
do
    outFile=${inFile%_toModify.png}.png
    echo "Input file : ${outFile}"
    echo "Cropping..."
    convert -background white -alpha remove ${inFile} -crop ${size}x0 ${outFile}

    for file in $(ls -1 ${outFile%.png}-?.png) ;
    do
	convert +repage ${file} _tmp ; mv _tmp ${file}
	identify ${file}
	echo "Trimming..."
	convert -trim ${file} ${file%.png}_trimmed.png
	convert +repage ${file%.png}_trimmed.png _tmp ; mv _tmp ${file%.png}_trimmed.png
	identify ${file%.png}_trimmed.png
	echo "Rotating..."
	convert -rotate 45 -background white -alpha remove ${file%.png}_trimmed.png _tmp_${file}
	convert +repage _tmp_${file} _tmp; mv _tmp _tmp_${file}
	identify _tmp_${file}
	echo "Cropping..."
	convert -background white -alpha remove _tmp_${file} -crop 1042x521-0-0 ${file%.png}_triang.png
	convert +repage ${file%.png}_triang.png _tmp ; mv _tmp ${file%.png}_triang.png
	identify ${file%.png}_triang.png
	rm -fr ${file} _tmp*png *_trimmed.png
    done
    rm -fr ${outFile} ${inFile}
done
