for f in $(find . -name "*.bed"); do
	#echo $f
	fname=$(echo $f | sed 's/.bed/.bed.bgz/')
	if [ ! -f fname ]; then
		echo $f
		cat $f | sort -k1,1 -k2,2n | awk '{split($4, a, ":"); $4 = a[1]; print}' | bgzip -c > $fname
	fi
done
