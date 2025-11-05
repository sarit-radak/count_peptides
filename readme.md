# Read Me

1. Make a project directory on the promethION
    1. `Cliff/Desktop/(name)` ?
2. Download [cyberduck](https://cyberduck.io/download/)
3. Connect to Joemetheus
    1. Requires VPN off-campus
    2. SFTP
    3. Server: 172.29.230.63
    4. Port: 22
    5. Username: cliff
    6. Password: newcomputer1
    
4. Bookmark project directory + data directory
    1. data directory: `/var/lib/minknow/data`
5. Download peptide counting code
    1. https://github.com/sarit-radak/count_peptides
6. Move directory to Joemetheus
7. Download [cursor](https://cursor.com)
8. Open directory in cursor
9. Run bash `count_peptides.sh`
10. Move this test file into `1-raw_data`
    
    [hpv_100k.fasta](Read%20Me/hpv_100k.fasta)
    
11. Uncomment the trim peptides step, comment out the count step
12. Run bash `count_peptides.sh`
    - Trim output (`logs/hpv_100k.log`)
        
        ```
        --- Summary ---
        Total reads processed: 100000
        Reads with both constants found (successful trims): 99670 (99.67%)
          Of successfully trimmed reads: 23923 (24.00%) are 24 bp
          Of successfully trimmed reads: 29774 (29.87%) are 27 bp
          Of successfully trimmed reads: 20813 (20.88%) are 30 bp
          Of successfully trimmed reads: 13865 (13.91%) are 33 bp
          Of successfully trimmed reads: 8646 (8.67%) are 36 bp
          Of successfully trimmed reads: 2649 (2.66%) are other lengths
        Reads with 5' missing (but 3' found): 26 (0.03%)
        Reads with 3' missing (but 5' found): 292 (0.29%)
        Reads with both missing: 12 (0.01%)
        Histogram (length -> count) written to: ./2-trimmed/hpv_100k_trimmed_length_hist.csv
        Examples (passed) saved to ./2-trimmed/hpv_100k_examples_passed.fasta
        Examples (missing 5') saved to ./2-trimmed/hpv_100k_examples_missing_5prime.fasta
        Examples (missing 3') saved to ./2-trimmed/hpv_100k_examples_missing_3prime.fasta
        Examples (missing both) saved to ./2-trimmed/hpv_100k_examples_missing_both.fasta
        
        Trimmed reads written to: ./2-trimmed/hpv_100k.fasta
        ```
        
13. Reverse commenting
    1. Leave `tolerate-errors` commented out of the count line
14. Move this list of peptide sequences into `reference`
    
    [hpv_l.fasta](Read%20Me/hpv_l.fasta)
    
15. Run bash `count_peptides.sh`
    - Count output
        
        ```
        --- Summary ---
        Sample: hpv_100k
        Total reads processed: 99670
        Reads matched (exact): 84749
        Reads matched (total): 84749 (85.03%)
        Reads unmatched: 14921 (14.97%)
        Peptide sequences: 571750
        Peptides identified (count > 0):   13233 (2.31%)
        Peptides not identified (count=0): 558517 (97.69%)
        ----------------
        Counts written to ./3-counts/hpv_100k.csv
        ```