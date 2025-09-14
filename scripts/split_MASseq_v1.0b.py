#!/home/zjia/miniconda3/bin/python
# This is the first release version, which enables a semi-automatic processing.
# This script demands a .json file caontaining the essential meta information.
# Current Version: v1.0 alpha

# Load the necessary libraries.
# Standard Python libraries:
import os
import sys
import time
import json
import getopt

# Third party packages:
import pysam
import regex

def getDatetime():
    return time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())

# ==================================== User Interface & Parameter Parsing ==================================== 
# Get the options provided by users in a dictionary.
usage = """This is the script to split MAS-ligated reads into "original transcripts" and validate split results. Valid results will output into a file with a ".BCassigned.tsv" extension. 
Invalid results will output into different files with different extensions, ".err.tsv" for CCS reads failed to split, ".deg.tsv" for split results without a "SigF" sequence in its 5' end, ".noBC.tsv" for split results without a identifiable 3' adapter barcode, ".noUMI.tsv" for split results without a detectable UMI pattern right after the "SigF" sequence.
This script will generate a json file containing some statistic information about its running process in the working directory by default.

General usage: 
  python split_MASseq_<version>.py [-p] [-h] [-m <meta_information_file>] [-v <valid_output_directory>] [-i <invalid_output_directory>] [<PATH>/]<file_name>.fastq

This script can either run within or without a standard project directory structure:
  -p    If this parameter is provided, this script will run in a "project mode" which needs a standard project directory structure;
  -m    The file name with or without the PATH to a .json files cantains necessary meta information, leave it NULL to find a "proj_meta.json" file in current WD;

The fillowing parameters is need when it is under a "standalone mode":
  -v    The directory to store the file containing valid reads, leave it NULL to output them in current WD;
  -i    The directory to store the file containing invalid reads for recall, leave it NULL to output them in current WD;

To view the usage information:
  -h    Print usage information and exit.
"""

optlist, args = getopt.getopt(sys.argv[1:], 'phm:v:i:f:')
optdict = dict(optlist)
projWD = os.getcwd()

if ("-h" in sys.argv[1:]) or (len(sys.argv)==1):
    sys.stderr.write(usage)
    sys.exit()


if len(args)==1:
    if args[0][-6:]==".fastq":
        fqf_name = os.path.basename(args[0][:-6])
    else:
        fqf_name = os.path.basename(args[0])
else:
    sys.stderr.write("This version could only process one .fastq file for each run :-(\nIf you want multi-threads processing, please refer to our manual and try 'ParaFly'.")
    sys.exit()

if '/' in fqf_name:
    fqf_name = fqf_name.split("/")

valid_dir = projWD
invalid_dir = projWD

if ("-v" in optdict.keys()) and optdict["-v"]:
    valid_dir = os.path.join(projWD, optdict["-v"])
if ("-i" in optdict.keys()) and optdict["-i"]:
    invalid_dir = os.path.join(projWD, optdict["-i"])


if "-p" in optdict.keys():
    print(f"[{getDatetime()}] Will run in a standard project directory structure, current WD: {projWD}")
    fq_file = f"{projWD}/fqsplit/{args[0]}"
    bca_file = f"{projWD}/valid/{fqf_name}.BCassigned.tsv"  # "BCA" for "BarCode Assigned".

    err_file = f"{projWD}/invalid/{fqf_name}.err.tsv"  # The sequence that failed to split.
    deg_file = f"{projWD}/invalid/{fqf_name}.deg.tsv"  # The 5' signature sequence was not intact.
    noBC_file = f"{projWD}/invalid/{fqf_name}.noBC.tsv"  # No 3'adapter barcode sequence was detected.
    noUMI_file = f"{projWD}/invalid/{fqf_name}.noUMI.tsv"  # No UMI pattern was detected.

    json_name = f"{projWD}/valid/{fqf_name}.stat.json"

else:
    print(f"[{getDatetime()}] Will run in a standalone mode, current WD: {projWD}")
    fq_file = os.path.join(projWD, args[0])
    bca_file = os.path.join(projWD, valid_dir, f"{fqf_name}.BCassigned.tsv")

    err_file = os.path.join(projWD, invalid_dir, f"{fqf_name}.err.tsv")
    deg_file = os.path.join(projWD, invalid_dir, f"{fqf_name}.deg.tsv")
    noBC_file = os.path.join(projWD, invalid_dir, f"{fqf_name}.noBC.tsv")
    noUMI_file = os.path.join(projWD, invalid_dir, f"{fqf_name}.noUMI.tsv")

    json_name = os.path.join(projWD, valid_dir, f"{fqf_name}.stat.json")


# ================================ Basic Information ====================================
# Load the meta information of the project. 
if ("-m" in optdict.keys()) and optdict["-m"]:
    meta_json = os.path.join(projWD, optdict["-m"])
else:
    meta_json = f"{projWD}/proj_meta.json"

print(f"[{getDatetime()}] The meta-information file will be: {meta_json}")

with open(meta_json) as metaf:
    meta_inf = json.load(metaf)

pattern_basic = {}
for i in meta_inf["UsedAdapter"]:
   pattern_basic[i] = "(?e)(" + meta_inf["AdapterBC"][i] + "){e<=2}"


# ================================ Defining File Handles ====================================
# Create handles of the output files according to the structure of the project.
# The structure of a standard splitting Project would be: 
# [GV/MII]_Project/
# ├── Orig_Data/
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
# └── scripts/
#     ├── project_meta.json
#     └── onestop.py
fq_bca = open(bca_file, "w")

fq_err = open(err_file, "w")  # The sequence that failed to split.
fq_deg = open(deg_file, "w")  # The 5' signature sequence was not intact.
fq_noBC = open(noBC_file, "w")  # No 3'adapter barcode sequence was detected.
fq_noUMI = open(noUMI_file, "w")  # No UMI pattern was detected.


# ================================= Defining Functions ====================================
# Loading meta information from a .json file.

# To get the complementary sequence of a given sequence. 
def seqComp(s):
	compDict = {'G':'C','C':'G','T':'A','A':'T','N':'N','-':'-'}
	return ''.join([compDict[x] for x in s])[::-1]


# Determine whether the read is positive or negative, to enable a unified workflow. 
# If positive, return its original sequence, else, return its complementary sequence.
def seqSigFWD(orig_seq, orig_qual):
    # Forward signature sequence, and the complementary sequence of reverse signature;
    fwd_det = bool(regex.search("(?e)(TCTACACGACGCTCTTCCGATCT){e<=2}", orig_seq)) and bool(regex.search("(?e)(CTCTGCGTTGATACCACTGCTTA){e<=2}", orig_seq))
    # Vice versa.
    bwd_det = bool(regex.search("(?e)(AGATCGGAAGAGCGTCGTGTAGA){e<=2}", orig_seq)) and bool(regex.search("(?e)(TAAGCAGTGGTATCAACGCAGAG){e<=2}", orig_seq))

    if fwd_det and (not bwd_det):
        return (orig_seq, orig_qual)

    elif (not fwd_det) and bwd_det:
        return (seqComp(orig_seq), orig_qual[::-1])
    
    else:
        return None
    

# Split the original reads with the complementary sequence of the signature sequence of the MAS primer (Reverse).
def splitPrim(fq_ID, fq_seq, fq_qual):
    ind_lis = [0]
    res_lis = []

    for hit in regex.finditer("(?e)(GTACTCTGCGTTGATACCACTGCTTA){e<=3}", fq_seq):
        ind_lis.extend(hit.span())

    for i in range(len(ind_lis)//2):
        res_lis.append((f"{fq_ID}|{i}", fq_seq[ind_lis[i*2]: ind_lis[i*2+1]], fq_qual[ind_lis[i*2]: ind_lis[i*2+1]]))

    return res_lis


# The function to check whether the 5' end of the split reads is intact or not.
def checkIntactSigF(sp_tup):
    sigf_hit = regex.search("(?e)(CTACACGACGCTCTTCCGATCT){e<=2}", sp_tup[1][:50])
    
    if sigf_hit:
        return (sp_tup[0], sp_tup[1][sigf_hit.span()[-1]:], sp_tup[-1][sigf_hit.span()[-1]:])
    else:
        return False
    

# The function to assign the 3'adapter BC for a given sequence. 
def adapterAssign(seq, pdic):
    for bc in pdic.keys():
        bc_hit = regex.search(pdic[bc], seq[-25:])
        if bc_hit:
            return (bc, bc_hit.span()[0]-25)

    return None
    

# ================================= Main ====================================
# Read the converted FastQ file using "pysam", and process by entry. 
# fq_file_handle = pysam.FastxFile(fq_file)
stat_dic = {
    "Split_failed": 0,  # CCS reads failed to split.
    "5end_deg": 0,  # Split reads without an intact 5' end.
    "No_BC": 0,  # Split reads without a detectable 3' adapter barcode.
    "No_UMI": 0,  # Split reads without a detectable UMI pattern. 
    "BC_assigned": 0  # Valid split reads.
}

for entry in pysam.FastxFile(fq_file):
    entry_ID = entry.name
    entry_seq = entry.sequence
    entry_qual = entry.quality

    lisP = seqSigFWD(entry_seq, entry_qual)

    if lisP:
        entry_seqP, entry_qualP = lisP
    else:
        fq_err.write(f"{entry_ID}|Error\t{entry_seq}\t{entry_qual}\n")
        stat_dic["Split_failed"] += 1
        continue

    sp_lis = splitPrim(entry_ID, entry_seqP, entry_qualP)
    for seq in sp_lis:
        ch_res = checkIntactSigF(seq)
        if ch_res:
            saa_res = adapterAssign(ch_res[1], pattern_basic)

            if saa_res:
                bca_ID = f"{ch_res[0]}|{saa_res[0]}"
                bca_seq = ch_res[1][:saa_res[-1]]
                bca_qual = ch_res[-1][:saa_res[-1]]

                matchUMI = regex.split("(^[ATCG]{8,12})(ATGGG){s<=1}", bca_seq, 1)

                if len(matchUMI)==4:
                    out_ID = f"{bca_ID}|{matchUMI[1]}"
                    out_seq = matchUMI[-1]
                    out_qual = bca_qual[-len(out_seq):]

                    fq_bca.write(f"{out_ID}\t{out_seq}\t{out_qual}\n")
                    stat_dic["BC_assigned"] += 1

                else:
                    fq_noUMI.write(f"{bca_ID}|noUMI\t{bca_seq}\t{bca_qual}\n")
                    stat_dic["No_UMI"] += 1
        
            else:
                fq_noBC.write(f"{ch_res[0]}|noBC\t{ch_res[1]}\t{ch_res[-1]}\n")
                stat_dic["No_BC"] += 1

        else:
            fq_deg.write(f"{seq[0]}|Degraded\t{seq[1]}\t{seq[-1]}\n")
            stat_dic["5end_deg"]

fq_bca.close()
fq_err.close()
fq_deg.close()
fq_noBC.close()
fq_noUMI.close()

print(f"[{getDatetime()}] All sequences sucessfully extracted :-)\nResult file: {bca_file}.")


# Dump statistic information into a .json file.
with open(json_name, "w") as jf:
    json.dump(stat_dic, jf, indent=4)

print(f"[{getDatetime()}] Json file: {json_name}.")