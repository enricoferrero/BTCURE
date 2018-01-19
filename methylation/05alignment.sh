#!/usr/bin/env bash

inDir="/GWD/bioinfo/projects/RD-PTS-EpinovaData-analysis/WGBS/02trimming"
outDir="/GWD/bioinfo/projects/RD-PTS-EpinovaData-analysis/WGBS/05alignment"
genomeDir="/GWD/bioinfo/projects/RD-PTS-EpinovaData-analysis/WGBS/00genome"
mkdir -p $outDir
rm -rf $outDir/*

cd $inDir
for sampleDir in Sample_*
do
	mkdir -p $outDir/$sampleDir
    cd $inDir/$sampleDir
    reads1=${sampleDir}-trimmed-pair1.fastq.gz
    reads2=${sampleDir}-trimmed-pair2.fastq.gz
    qsub -cwd -V -q epigen -pe epigen 8 -N ${sampleDir}"_05" -hold_jid "*_04" -o $outDir/bismark.log.txt -e $outDir/bismark.log.txt -b yes \
        bismark --bowtie2 -p 8 --temp_dir $outDir/$sampleDir --output_dir $outDir/$sampleDir $genomeDir -1 $reads1 -2 $reads2
done
