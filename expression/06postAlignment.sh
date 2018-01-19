#!/usr/bin/env bash

inputDir=/home/ef884766/Projects/RASynovial/data/analysis/RNAseq/05alignment
outputDir=/home/ef884766/Projects/RASynovial/data/analysis/RNAseq/09postAlignment

cd $inputDir
for i in Sample_RA-*/Aligned.out.sam
do
	# grab sample name
	sampleName=$(dirname $i)
	# remove unmapped reads (0x4 = 4)
	# remove multiple mappings/secondary alignments (0x100 = 256)
	# 4 + 256 = 260 (-F 260)
	# remove reads with mapping score < 10 (-q 10)
	# remove mitochondrial sequences (awk '($1 ~ /^@/) || ($3 != "MT") { print $0 }')
	# convert SAM to BAM (-bS)
	samtools view -@ 8 -S -h -F 260 -q 10 $i | awk '($1 ~ /^@/) || ($3 != "MT") { print $0 }' | samtools view -@ 8 -bS -o $outputDir/${sampleName}.tmp.bam - |& tee -a $outputDir/samtools.log.txt
	# sort reads
	samtools sort -@ 8 -m 4G $outputDir/${sampleName}.tmp.bam $outputDir/$sampleName |& tee -a $outputDir/samtools.log.txt
	# create index
	samtools index $outputDir/${sampleName}.bam |& tee -a $outputDir/samtools.log.txt
	# clean up
	rm $outputDir/${sampleName}.tmp.bam
done
