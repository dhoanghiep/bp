bedtools intersect -wa -wb -a ../bed_file/introns.bed -b ../bed_file/mercer.bed | awk '{printf "%s\t%i\t%i\t%s\t%i\t%s\t%i\n", $1,$2,$3,$4,$8,$6,$11}'  > ../bed_file/introns_with_bp_mercer.bed; 

bedtools intersect -wa -wb -a ../bed_file/introns.bed -b ../bed_file/pineda.bed | awk '{printf "%s\t%i\t%i\t%s\t%i\t%s\t%i\n", $1,$2,$3,$4,$8,$6,$11}'  > ../bed_file/introns_with_bp_pineda.bed