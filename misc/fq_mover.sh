cat controls_metadata.tsv | grep -v Lane | while read lyne; do 
	barcode=$(echo $lyne | cut -f 2 -d " " );
	otherName=$(echo $lyne | cut -f 1 -d " " );
	mkdir -p /proj/cdjones_lab/csoeder/EJgrepper/FASTQ/controls/"$otherName"/ ; 
	echo "cat /proj/cdjones_lab/csoeder/EJgrepper/FASTQ/MOLNG-2416/n_[1-4]_1_"$barcode".fastq > /proj/cdjones_lab/csoeder/EJgrepper/FASTQ/controls/"$otherName"/"$otherName".R1.fq"
	echo "cat /proj/cdjones_lab/csoeder/EJgrepper/FASTQ/MOLNG-2416/n_[1-4]_2_"$barcode".fastq > /proj/cdjones_lab/csoeder/EJgrepper/FASTQ/controls/"$otherName"/"$otherName".R2.fq"
done | awk '{print "|"$0"|"}' | tr '|' '"' | awk '{print "sbatch -n 1 -t 1:00:00 --mem=8G --wrap="$0}'
