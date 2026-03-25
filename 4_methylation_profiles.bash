#!/usr/bin/env bash
###This function takes a BAM file of ONT reads containing methylation tags MM/ML, processes it, and outputs a trimmed BAM ready to be re-aligned.
modbam_processing() {
    local out=$out

	# sort bam by read name
	samtools sort \ 
	 -@ 1 \
	 -o tmp_fastq/${out}.sort.bam \
	 -n \
	 ${file}

	# convert to fastq
	samtools fastq \
	 -@ 1 \
	 -T MM,ML \
	 tmp_fastq/${out}.sort.bam > tmp_fastq/${out}.fastq

	# run porechop to trim adapters
	/projects/novabreed/share/gmagris/software/Porechop/./porechop-runner.py \
	 -i tmp_fastq/${out}.fastq \
	 -t 1 \
	 --no_split \
	 --extra_end_trim 10 \
	 --check_reads 500 \
	 -o tmp_fastq/${out}.trim.fastq

#Convert trimmed FASTQ back to BAM
	samtools import \
	 -@ 1 \
	 tmp_fastq/${out}.trim.fastq \
	 -T MM,ML \
	 -o tmp_fastq/${out}.trim.bam 

	rm -rf tmp_fastq/${out}.fastq tmp_fastq/${out}.trim.fastq

#Repair methylation tags with Modkit (MM/ML tags encode methylated positions relative to the original sequence. 
#-->modkit repair recomputes MM/ML so that they match the trimmed reads.
	../dist_modkit_v0.5.1_8fa79e3/modkit repair \
	-t 1 \
	-d tmp_fastq/${out}.sort.bam \
	-a tmp_fastq/${out}.trim.bam \
	-o tmp_aln/${out}.trim.bam \
	--log-filepath ./${out}_repair.log
	
	rm -rf tmp_fastq/${out}.sort.bam tmp_fastq/${out}.trim.bam
}

mkdir -p tmp_aln tmp_fastq tmp_aln_dorado
max_jobs=20

sample_list=("Boni" "T116" "T206" "T288" "T33" "T345" "T379" "T87" "Dtt2" "T132" "T224" "T297" "T330" "T365" "T404" "T99" "T177" "T283" "T304" "T336" "T366" "T412")
# Quncho  

for sample in "${sample_list[@]}"
do
	# launch background jobs, limited to $max_jobs
	for file in ../ont_bam_files/${sample}/*.bam
	do
		out=$(basename ${file} .bam)
		# control concurrency
		while (( $(jobs -rp | wc -l) >= max_jobs )); do
			sleep 1
		done

		# launch job in background
		modbam_processing $out &
