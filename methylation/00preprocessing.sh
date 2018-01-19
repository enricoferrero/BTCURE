#!/usr/bin/env bash

inDir="/GWD/bioinfo/projects/RD-PTS-EpinovaData-raw/WGBS"
outDir="/GWD/bioinfo/projects/RD-PTS-EpinovaData-analysis/WGBS/00preprocessing"
mkdir -p $outDir
rm -rf $outDir/*

cd $inDir
for dir in Project_Methyl_Batch*_Run1_*
do
    cd $dir
    md5sum -c RunStat/*_md5.txt |& tee -a $outDir/md5sum.log.txt &
    cd $inDir
    
    for sample in $dir/Sample_*
    do
        sample=$(basename $sample)
        mkdir -p $outDir/$sample
        
        for file in $dir/$sample/*.fastq.gz
        do
            file=$(basename $file)
            ln -s -f $inDir/$dir/$sample/$file $outDir/$sample/${file/L0/L1}
        done

    done

done

cd $inDir
for dir in Project_Methyl_Batch*_Run2_*
do
    cd $dir
    md5sum -c RunStat/*_md5.txt |& tee -a $outDir/md5sum.log.txt &
    cd $inDir
    
    for sample in $dir/Sample_*
    do
        sample=$(basename $sample)
        mkdir -p $outDir/$sample
        
        for file in $dir/$sample/*.fastq.gz
        do
            file=$(basename $file)
            ln -s -f $inDir/$dir/$sample/$file $outDir/$sample/${file/L0/L2}
        done

    done

done

wait

cd $outDir
for sampleDir in Sample_*
do
    echo "unpigz -c $sampleDir/*_R1_*.fastq.gz | pigz --fast > $sampleDir/${sampleDir}_R1.fastq.gz" | \
        qsub -cwd -V -q epigen -pe epigen 8 -N ${sampleDir}"_R1_00" -o gzip.log.txt -e gzip.log.txt
    echo "unpigz -c $sampleDir/*_R2_*.fastq.gz | pigz --fast > $sampleDir/${sampleDir}_R2.fastq.gz" | \
        qsub -cwd -V -q epigen -pe epigen 8 -N ${sampleDir}"_R2_00" -o gzip.log.txt -e gzip.log.txt
done
