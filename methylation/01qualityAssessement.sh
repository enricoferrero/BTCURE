#!/usr/bin/env bash

inDir="/GWD/bioinfo/projects/RD-PTS-EpinovaData-analysis/WGBS/00preprocessing"
outDir="/GWD/bioinfo/projects/RD-PTS-EpinovaData-analysis/WGBS/01qualityAssessment"
mkdir -p $outDir
rm -rf $outDir/*

cd $inDir
for sampleDir in Sample_*
do
	qsub -cwd -V -q epigen -pe epigen 2 -N ${sampleDir}"_01" -hold_jid "*_00" -o $outDir/fastqc.log.txt -e $outDir/fastqc.log.txt -b yes \
        fastqc --outdir=$outDir --threads 2 $sampleDir/${sampleDir}_R[12].fastq.gz
done
