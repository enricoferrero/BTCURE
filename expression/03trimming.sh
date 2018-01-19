#!/usr/bin/env bash

rawDataDir=/GWD/bioinfo/projects/RD-PTS-EpinovaData-raw/RNAseq
resultsDir=/GWD/bioinfo/projects/RD-PTS-EpinovaData-analysis/RNAseq/03trimming

cd $rawDataDir
for sampleDir in Sample_RA-*
do
	mkdir -p $resultsDir/$sampleDir
	cd $rawDataDir/$sampleDir
	# creates arrays for paired-end reads
	reads1=(*_R1_*.fastq.gz)
	reads2=(*_R2_*.fastq.gz)
	# loop over array(s) indices
	for (( i = 0; i < ${#reads1[@]}; i++ ))
	do
		# scythe: adapter removal
		cd $rawDataDir/$sampleDir
		scythe -a $resultsDir/IlluminaTruSeqAdapters.fa -o $resultsDir/$sampleDir/${reads1[$i]%.fastq.gz}.tmp.fastq ${reads1[$i]} -M 20 |& tee -a $resultsDir/scythe.log.txt &
		scythe -a $resultsDir/IlluminaTruSeqAdapters.fa -o $resultsDir/$sampleDir/${reads2[$i]%.fastq.gz}.tmp.fastq ${reads2[$i]} -M 20 |& tee -a $resultsDir/scythe.log.txt &
		wait
		# sickle: quality trimming
		cd $resultsDir/$sampleDir
		sickle pe -t sanger -g -f ${reads1[$i]%.fastq.gz}.tmp.fastq -r ${reads2[$i]%.fastq.gz}.tmp.fastq -o ${reads1[$i]} -p ${reads2[$i]} -s ${sampleDir}.singles.fastq.gz -q 20 -l 20 |& tee -a $resultsDir/sickle.log.txt &
	done
	rm *.tmp.fastq
done
