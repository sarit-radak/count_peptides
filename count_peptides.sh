#!/bin/bash
mkdir -p logs
mkdir -p 1-raw_data
mkdir -p 2-trimmed
mkdir -p 3-counts

libraries=($(find 1-raw_data -mindepth 1 -maxdepth 1 -type f -exec bash -c 'for f; do basename "${f%.fasta}"; done' _ {} +))

for library in "${libraries[@]}"; do
    log_file="logs/${library}.log"
    {

    ### trim reads
    #python3 -u ./scripts/trim.py -r ./1-raw_data/$library.fasta -c ./reference/constant.fasta -o ./2-trimmed/$library.fasta --troubleshoot

    ### count peptides
    python3 -u ./scripts/count.py -r ./2-trimmed/$library.fasta -p ./reference/hpv_l.fasta -o ./3-counts/$library.csv #--tolerate-errors --threads 1 --batch-size 2000 --min-shared-kmers 5 #-t 16 

    } >"$log_file" 2>&1 &

done

# wait for all backgrounded jobs to finish
wait