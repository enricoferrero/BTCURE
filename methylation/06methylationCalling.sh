#!/usr/bin/env bash

inDir="/GWD/bioinfo/projects/RD-PTS-EpinovaData-analysis/WGBS/05alignment"
outDir="/GWD/bioinfo/projects/RD-PTS-EpinovaData-analysis/WGBS/06methylationCalling"
mkdir -p $outDir
rm -rf $outDir/*

cd $inDir
for sample in Sample_*/*.bam
do
	sampleName=$(dirname $sample)
    mkdir -p $outDir/$sampleName
    qsub -cwd -V -q epigen -pe epigen 8 -N ${sampleName}"_06" -hold_jid "*_05" -o $outDir/bismark.log.txt -e $outDir/bismark.log.txt -b yes \
        bismark_methylation_extractor --multicore 8 --paired-end --comprehensive --merge_non_CpG --no_header --gzip --bedGraph --zero_based --ample_memory --output $outDir/$sampleName $sample
done
