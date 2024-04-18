#!/bin/bash

FASTQ_FILE="/home/mariopasc/C++/bioinformatics-algorithms/fastq_files/ERR103404_1.fastq.gz"

./build/main kmerfreq "$FASTQ_FILE" > workflow_scripts/temp_output.txt
python ./Assembly/De_Brujin_Graphs/visualize_kmerFrequency.py < workflow_scripts/temp_output.txt
