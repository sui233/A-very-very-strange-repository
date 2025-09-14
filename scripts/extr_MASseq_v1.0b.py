# Author: JIA Zheng
# Created on 2025/07/26
# This is the script to initialize a MAS-PAIso-seq2 data splitting project.
# This script will extract the sequence, quality, ID and pass number of every reads into a FastQ File.
# Current version: 1.0-beta

# Load the necessary libraries.
# Standard Python libraries:
import os
import sys
import getopt
import time
import json

# Third party packages:
import pysam

def getDatetime():
    return time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())

# ==================================== User Interface & Parameter Parsing ==================================== 
# Get the options provided by users in a dictionary.
usage = """This is a script to extract sequences with its ID, quality information and pass number from a .bam file containing MAS-PAIso-seq(2) sequencing information to a .fastq file.
This script will generate a json file containing some statistic information about its running process in the working directory by default.

General usage: 
  python extr_MASseq_<version>.py [-p] [-h] [-o <output_directory>] [-f <output_filename>] [<PATH>/]<file_name>.bam

This script can either run within or without a standard project directory structure:
  -p    If this parameter is provided, this script will run in a "project mode" which needs a standard project directory structure;

The fillowing parameters is needed when it is under a "standalone mode":
  -o    The PATH to which the output files will be created, leave it NULL to output them in current WD;
  -f    The name for output files, leave it null to use the default name.

To view the usage information:
  -h    Print usage information and exit.
"""

optlist, args = getopt.getopt(sys.argv[1:], 'hpo:f:')
optdict = dict(optlist)
projWD = os.getcwd()

if ("-h" in sys.argv[1:]) or (len(sys.argv)==1):
    sys.stderr.write(usage)
    sys.exit()

# Load meta information the project.
# Configuration for the "project" mode.
if "-p" in optdict.keys():
    print(f"[{getDatetime()}] Will run in a 'project' mode, current WD: {projWD}")
    
    for file_name in os.listdir("hifi_reads/"):
        if file_name[-4:]==".bam":
            bamf_name = f"{projWD}/hifi_reads/{file_name}"
            break

    print(f"[{getDatetime()}] BAM file {file_name} will be processed.")

    gt3_name = f"{projWD}/css.fastq"
    lt3_name = f"{projWD}/discard/css_pass_lt3.fastq"
    err_sam_name = f"{projWD}/recall/potential_err.sam"

    json_name = f"{projWD}/passnum_stat.json"
            
else:
    print(f"[{getDatetime()}] Will run in a 'standard' mode, current WD: {projWD}")
    if len(args)==1:
        bamf_name = args[0]
    else:
        sys.stderr.write("The name of the BAM file should be provided in a standalone mode.")
        sys.exit()

    otp_dir = projWD
    otp_fname = "css"

    if ("-o" in optdict.keys()) and optdict["-o"]:
        otp_dir = os.path.join(projWD, optdict["-o"])

    if ("-f" in optdict.keys()) and optdict["-f"]:
        otp_fname = optdict["-f"]

    otp_name = os.path.join(otp_dir, otp_fname)

    gt3_name = f"{otp_name}.fastq"
    lt3_name = f"{otp_name}_pass_lt3.fastq"
    err_sam_name = f"{otp_name}_no_passnum.sam"

    json_name = f"{otp_name}_passnum_stat.json"

print(f"[{getDatetime()}] {bamf_name} will be processed.")

# ==================================== Defining Functions ==================================== 
# Loading meta information from a .json file.

# The Function to determine the index of the "pass number" information.
def getPassIndex(bamf):
    """
    The function to get the colunm index of 'pass number' from the first line.
    
    Arg:
      bamf (file handle): the handle of the .bam file, which is provided by the 'AlignmentFile' file from 

    Return: return_description
    """
    
    with pysam.AlignmentFile(bamf, "rb", check_sq=False) as bf:
        frec = next(bf)

    rec_dic = frec.to_dict()
        
    for i in range(len(rec_dic["tags"])):
        if rec_dic["tags"][i][:5] == "np:i:":
            return i


# ==================================== Main loop ====================================
pn_ind = getPassIndex(bamf_name)
print(f"[{getDatetime()}] The default colunm index of 'pass number' is set on: {pn_ind}.")

gt3_fq = open(gt3_name, "w")
lt3_fq = open(lt3_name, "w")
err_sam = open(err_sam_name, "w")

stat_dic = {
    "PNgt3": 0, 
    "PNlt3": 0, 
    "otherPNcol": 0,
    "noPN": 0
}

for query in pysam.AlignmentFile(bamf_name, "rb", check_sq=False):
    samq_dict = query.to_dict()

    if samq_dict['tags'][pn_ind][:5] == "np:i:":
        pass_num = int(samq_dict['tags'][pn_ind][5:])

    else:
        no_pn = True
        for i in range(len(samq_dict['tags'])):
            if samq_dict['tags'][i][:5] == "np:i:":
                pass_num = int(samq_dict['tags'][i][5:])
                no_pn = False

                print(f"[{getDatetime()}] Sequence entry '{samq_dict['name']}' has its pass number information in column {i}.")
                stat_dic["otherPNcol"] += 1
                break

        if no_pn:
            err_sam.write(f"{query.to_string()}\n")

            print(f"[{getDatetime()}] Sequence entry '{samq_dict['name']}' doesn't have pass number information.")
            stat_dic["noPN"] += 1
            continue

    if pass_num>=3:
        gt3_fq.write(f"@{samq_dict['name']}|{pass_num}\n{samq_dict['seq']}\n+\n{samq_dict['qual']}\n")
        stat_dic["PNgt3"] += 1

    else:
        lt3_fq.write(f"@{samq_dict['name']}|{pass_num}\n{samq_dict['seq']}\n+\n{samq_dict['qual']}\n")
        stat_dic["PNlt3"] += 1


gt3_fq.close()
lt3_fq.close()
err_sam.close()

print(f"[{getDatetime()}] All sequences sucessfully extracted :-)\nResult file: {gt3_name}.")


# Dump statistic information into a .json file.
with open(json_name, "w") as jf:
    json.dump(stat_dic, jf, indent=4)

print(f"[{getDatetime()}] Json file: {json_name}.")