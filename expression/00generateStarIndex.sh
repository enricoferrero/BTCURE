#!/usr/bin/env bash
genomeDir=/GWD/bioinfo/projects/RD-PTS-EpinovaData-analysis/RNAseq/00genome

STAR --runThreadN 8 --runMode genomeGenerate --genomeDir $genomeDir --genomeFastaFiles $genomeDir/Homo_sapiens.GRCh38.dna.primary_assembly.fa --sjdbGTFfile $genomeDir/Homo_sapiens.GRCh38.78.gtf --sjdbOverhang 99
