cat ../hg19/gencode.v19.annotation.gtf | python3 extract_intron.py | sort -u -k1,1 -k2,3n > ../bed_file/introns.bed
