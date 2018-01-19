#!usr/bin/env bash

source_dir="/home/al580162/projects/BTCURE/Secondary_analyses"

data_dir="/GWD/bioinfo/projects/RD-PTS-EpinovaData-analysis/RNAseq/00genome"
reference_fa=$data_dir"/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
reference_gtf=$data_dir"/Homo_sapiens.GRCh38.78.gtf"

output_dir=$source_dir"/output/05_Cellular_enrichment/rsem_reference"

rsem-prepare-reference --gtf $reference_gtf $reference_fa $output_dir
