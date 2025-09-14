# 2025/09/03
# Author: JIA Zheng
# This is the script to separate the result of the splitting workflow. 
# The output files are named after their sample names or their barcodes.
# The current version of this script only works in a "project" mode. 
# This script demands a .json file caontaining the essential meta information. 
# nohup python convert_tsv2fqgz.py &
# Current version: 1.0-beta


import os
import time

import gzip 
import json 

# ================================= Defining Functions ====================================
# Get current date time.
def getDatetime():
    return time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())


# ================================ Basic Information Loading ====================================
# Loading meta information from a .json file.
with open('proj_meta.json') as metaf:
    meta_inf = json.load(metaf)

gz_handle_dic = {}
rnum_stat = {}

# If there isn't a "Adapter2Sample" key in the meta information file, the output files will be named after their barcodes. 
samp_name_dic = {}
if "Adapter2Sample" in meta_inf.keys():
    for bc in meta_inf["UsedAdapter"]:
        samp_name_dic[bc] = meta_inf["Adapter2Sample"][bc]
    else:
        samp_name_dic[bc] = bc

for bc in meta_inf["UsedAdapter"]:
    gz_handle_dic[bc] = gzip.open(f'split_result/{samp_name_dic[bc]}.fastq.gz', 'wt', compresslevel=6)
    rnum_stat[bc] = 0

conv_files = [f for f in os.listdir("valid/") if f[-4:]==".tsv"]


# ================================= Main loop ====================================
print(f"[{getDatetime()}] Timer started.")
curr_time = time.time()
start_time = curr_time

for tsvf in conv_files:
    print(f"Converting .tsv file {tsvf}...", end='\t')
    for line in open(f"valid/{tsvf}"):
        entryLis = line.strip().split()
        gz_handle_dic[entryLis[0].split("|")[3]].write(f"@{entryLis[0]}\n{entryLis[1]}\n+\n{entryLis[-1]}\n")
        rnum_stat[entryLis[0].split("|")[3]] += 1
    print(f"Time cost: {time.time() - curr_time}")
    curr_time = time.time()

print(f"[{getDatetime()}] Totally cost: {curr_time - start_time}s.")

for bc in meta_inf["UsedAdapter"]:
    gz_handle_dic[bc].close()


# Dump statistic information into a .json file.
with open("sample_reads.stat") as f:
    for bc in meta_inf["UsedAdapter"]:
        f.write(f"{bc}\t{meta_inf["Adapter2Sample"][bc]}\t{rnum_stat[bc]}\n")