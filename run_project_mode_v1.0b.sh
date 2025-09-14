#!/usr/bin/bash
# This is the bash script to run scripts in a "project" mode, which will automatically finish the split and recall process of the MAS-PAIso-seq(2) data.
# This script needs a "proj_meta.json" file that contains basic meta-information of this project to run properly. 
# To prepare a "proj_meta.json" file, please refer to our "user manual".

# Initiating a standard project structure for the project. 
# If there is a 
# The structure of a standard project would be: 
# [PROJ_name]_Project/
# ├── hifi_reads/
# │   ├── <filename>.bam
# │   ├── <filename>.bam.pbi
# │   └── <filename>.bam.xml
# ├── css.fastq
# ├── recall/
# ├── discarded/
# ├── fqsplit/
# │   ├── css.fastq-partxxx
# │   └── ...
# ├── valid/
# ├── invalid/
# ├── run_project_mode_v1.0b.sh
# ├── proj_meta.json
# └── scripts/
#     ├── extr_MASseq_v1.0b.py
#     ├── fastq-splitter.pl
#     ├── split_MASseq_v1.0b.py
#     ├── recall_MASseq_v1.0b.py
#     ├── false_split_detect_v1.0b.py
#     └── convert_tsv2fqgz_v1.0b.py
mkdir fqsplit
mkdir valid
mkdir invalid
mkdir recall
mkdir recall/false_split
mkdir discard
mkdir split_result

mkdir log_files
mkdir log_files/fqsplit_log

# Extracting CCS reads and their pass number from the .bam file.
python -u scripts/extr_MASseq_v1.0a.py -p > log_files/seq_extract.log

scripts/fastq-splitter.pl --n-parts 100 css.fastq
mv css.part-*.fastq fqsplit/
rm css.fastq

# Splitting CCS reads in 50 threads.
for i in {01..100}; do echo "python -u scripts/split_MASseq_v1.0a.py -p css.part-${i}.fastq > log_files/fqsplit_log/part-${i}_split.log"; done > mas_split.sh 
ParaFly -c mas_split.sh -CPU 50

python -u scripts/recall_MASseq_v1.0a.py -p 

# Finding potential false splits, which would not take a long time.
python -u scripts/false_split_detect_v1.0a.py -p

cat valid/*.tsv | awk -F '|' '{print $4}' | sort | uniq -c > sample_reads.txt

python -u scripts/convert_tsv2fqgz_v1.0a.py > log_files/compress_by_sample.log
