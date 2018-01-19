#!/usr/bin/env bash
genomeDir=/GWD/bioinfo/projects/RD-PTS-EpinovaData-analysis/WGBS/00genome/lambda
mkdir -p $genomeDir

bismark_genome_preparation --verbose --bowtie2 $genomeDir |& tee -a $genomeDir/bismark_genome_preparation.log.txt
