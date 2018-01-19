#!/usr/bin/env bash

inDir="/GWD/bioinfo/projects/RD-PTS-EpinovaData-analysis/WGBS/00preprocessing"
outDir="/GWD/bioinfo/projects/RD-PTS-EpinovaData-analysis/WGBS/02trimming"
mkdir -p $outDir
rm -rf $outDir/*

cd $inDir
for sampleDir in Sample_*
do
	mkdir -p $outDir/$sampleDir
	cd $inDir/$sampleDir
	reads1=${sampleDir}_R1.fastq.gz
	reads2=${sampleDir}_R2.fastq.gz
	qsub -cwd -V -q epigen -pe epigen 8 -N ${sampleDir}"_02" -hold_jid "*_01" -o $outDir/skewer.log.txt -e $outDir/skewer.log.txt -b yes \
        skewer --threads 8 --output $outDir/$sampleDir/${reads1%_R1.fastq.gz} --end-quality 20 --mean-quality 20 --min 20 --compress $reads1 $reads2
done
