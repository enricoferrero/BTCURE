#!usr/bin/env bash

source_dir="/home/al580162/projects/BTCURE/Secondary_analyses"
txome_dir=$source_dir/"output/05_Cellular_enrichment/rsem_reference"
output_dir=$source_dir/"output/05_Cellular_enrichment/salmon_alignment"

salmon index -t $txome_dir/"GRCh38.transcripts.fa" -i $output_dir/"GRCh38.transcripts_index" --type quasi -k 31
