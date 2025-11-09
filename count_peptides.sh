#!/bin/bash
mkdir -p logs
mkdir -p 1-raw_data
mkdir -p 2-trimmed
mkdir -p 3-counts


# copy fastq
#cp -r /var/lib/minknow/data/20251106_Jonah_and_Markus_library/no_sample_id/20251106_1531_MN49129_FBC15226_34bd3843/fastq_pass/ ./fastq_pass

# trim reference
#python3 -u ./scripts/trim.py -r ./reference/Armstrong_CO_peptides.fasta -c ./reference/constant.fasta -o ./reference/Armstrong_CO_peptides_trimmed.fasta --troubleshoot



libraries=($(find 1-raw_data -mindepth 1 -maxdepth 1 -type f -exec bash -c 'for f; do basename "${f%.fasta}"; done' _ {} +))

for library in "${libraries[@]}"; do
    log_file="logs/${library}.log"
    {
    echo $log_file
    ### trim reads
    #python3 -u ./scripts/trim.py -r ./1-raw_data/$library.fasta -c ./reference/constant.fasta -o ./2-trimmed/$library.fasta --troubleshoot

    ### count peptides
    #python3 -u ./scripts/count.py -r ./2-trimmed/$library.fasta -p ./reference/Armstrong_CO_peptides_trimmed.fasta -o ./3-counts/$library.csv #--tolerate-errors --threads 1 --batch-size 2000 --min-shared-kmers 5 #-t 16 


    } >"$log_file" 2>&1 &

done

# wait for all backgrounded jobs to finish
wait

python3 -u scripts/collect_frequencies.py