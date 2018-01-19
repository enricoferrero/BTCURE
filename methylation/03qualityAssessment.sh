#!/usr/bin/env bash

inDir="/GWD/bioinfo/projects/RD-PTS-EpinovaData-analysis/WGBS/02trimming"
outDir="/GWD/bioinfo/projects/RD-PTS-EpinovaData-analysis/WGBS/03qualityAssessment"
mkdir -p $outDir
rm -rf $outDir/*

cd $inDir
for sampleDir in Sample_*
do
	qsub -cwd -V -q epigen -pe epigen 2 -N ${sampleDir}"_03" -hold_jid "*_02" -o $outDir/fastqc.log.txt -e $outDir/fastqc.log.txt -b yes \
        fastqc --outdir=$outDir --threads 2 $sampleDir/*.fastq.gz
done
