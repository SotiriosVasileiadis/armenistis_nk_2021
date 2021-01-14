for i in `cat mylist`
	do
	TrimmomaticPE -threads 4 -phred33 ../../2_demultiplex_data/bac/NP2_out/lOTUs_moth_basis/${i}_lib_R1.fastq ../../2_demultiplex_data/bac/NP2_out/lOTUs_moth_basis/${i}_lib_R2.fastq trimmed/${i}_tr_R1.fastq trimmed/${i}_tr_R1_unpaired.fastq trimmed/${i}_tr_R2.fastq trimmed/${i}_tr_R2_unpaired.fastq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:100
	flash -m 10 -M 250 -x 0.1 -p 33 -d combined_merged -o ${i} --compress-prog=gzip --suffix=gz -t 4 trimmed/${i}_tr_R1.fastq trimmed/${i}_tr_R2.fastq
	done
