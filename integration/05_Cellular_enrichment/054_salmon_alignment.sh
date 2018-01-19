#!/usr/bin/env bash

source_dir="/home/al580162/projects/BTCURE/Secondary_analyses"
rawdata_dir="/GWD/bioinfo/projects/RD-PTS-EpinovaData-analysis/RNAseq/03trimming"
salmonReference=$source_dir/"output/05_Cellular_enrichment/salmon_reference/GRCh38.transcripts_index"
output_dir=$source_dir/"output/05_Cellular_enrichment/salmon_alignment"


cd $rawdata_dir
for sample in Sample_RA-*
do
	mkdir -p $outDir/$sampleDir
	read1=$(ls $rawdata_dir/$sample/RA-*_R1_*.fastq.gz)
	read2=$(ls $rawdata_dir/$sample/RA-*_R2_*.fastq.gz)
	
	#STAR --runThreadN 8 --genomeDir $genomeDir --readFilesIn $read1 $read2 --readFilesCommand zcat --outFileNamePrefix $resultsDir/$sampleDir/ |& tee -a $resultsDir/STAR.log.txt
	salmon quant -i $salmonReference -l A -1 $read1 -2 $read2 -o $output_dir/$sample
done
