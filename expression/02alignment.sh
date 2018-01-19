#!/usr/bin/env bash

rawDataDir=/GWD/bioinfo/projects/RD-PTS-EpinovaData-raw/RNAseq
genomeDir=/GWD/bioinfo/projects/RD-PTS-EpinovaData-analysis/RNAseq/00genome
resultsDir=/GWD/bioinfo/projects/RD-PTS-EpinovaData-analysis/RNAseq/02alignment

cd $rawDataDir
for sampleDir in Sample_RA-*
do
	mkdir $resultsDir/$sampleDir
	read1=$(ls -m $rawDataDir/$sampleDir/RA-*_R1_*.fastq.gz | tr -d '\n' | tr -d ' ')
	read2=$(ls -m $rawDataDir/$sampleDir/RA-*_R2_*.fastq.gz | tr -d '\n' | tr -d ' ')
	
	STAR --runThreadN 8 --genomeDir $genomeDir --readFilesIn $read1 $read2 --readFilesCommand zcat --outFileNamePrefix $resultsDir/$sampleDir/
done
