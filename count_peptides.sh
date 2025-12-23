#!/bin/bash

# make directories
mkdir -p logs
mkdir -p 1-input_fasta
mkdir -p 2-trimmed
mkdir -p 3-counts



# trim reference
#python3 -u ./scripts/trim.py -r ./reference/SIV_Pooled_Libraries.fasta -c ./reference/constant.fasta -o ./reference/SIV_Pooled_Libraries_trimmed.fasta

# rename peptides that appear multiple times
#python3 -u ./scripts/rename_duplicates.py -i ./reference/SIV_Pooled_Libraries.fasta -o ./reference/SIV_Pooled_Libraries_renumbered.fasta



# copy fastq_pass directory
#cp -r /var/lib/minknow/data/(path)/fastq_pass/ ./fastq_pass

# rename data folders (barcode01->Naive)
#python3 -u scripts/rename_folders.py



# pull library names from fastq folders OR
#libraries=($(find ./fastq_pass -mindepth 1 -maxdepth 1 -type d -exec bash -c 'for f; do basename "${f%.fasta}"; done' _ {} +))

# pull library names from fasta files
libraries=($(find 1-input_fasta -mindepth 1 -maxdepth 1 -type f -exec bash -c 'for f; do basename "${f%.fasta}"; done' _ {} +))

for library in "${libraries[@]}"; do
    log_file="logs/${library}.log"
    {
    echo "Processing: "$library

    # collect fastq files into single fasta
    #python3 -u ./scripts/collect_to_fasta.py -l $library -i ./fastq_pass/$library -o ./1-input_fasta
    
    # trim reads
    #python3 -u ./scripts/trim.py -r ./1-input_fasta/$library.fasta -c ./reference/constant.fasta -o ./2-trimmed/$library.fasta

    # count peptides
    #python3 -u ./scripts/count.py -r ./2-trimmed/$library.fasta -p ./reference/test_ref_renumbered.fasta -o ./3-counts/$library.csv

    } >"$log_file" 2>&1 &

done

# wait for all backgrounded jobs to finish
wait

# collect summary data for trim and count steps
#python3 -u scripts/collect_summaries.py

# collect frequencies across libraries
#python3 -u scripts/collect_frequencies.py
