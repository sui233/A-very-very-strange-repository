# 2025/08/13
# Author: JIA Zheng
# This is a script to identify false splits in the split results (entries) to provide a basis for downstream recalling.
# This script will provide a candidate list of false splits.
# This script is optional for the workflow and is suggested to run under the "project" mode.
# Current version: 1.0-beta

# Load the necessary libraries.
# Standard Python libraries:
import os
import sys
import getopt
import time
import json

# Third party packages:
import regex

def getDatetime():
    return time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())

# ==================================== User Interface & Parameter Parsing ==================================== 
# Get the options provided by users in a dictionary.
usage = """This is a script to detect potential false splits conducted by CCS reads splitting.
When this script run in a "standalone" mode, its input files should be prepared manually. 
This script will generate a json file containing some statistic information about its running process in the working directory by default.

General usage: 
  python split_MASseq_<version>.py [-p] [-h] [-m <meta_information_file>] [-c <candidate_output_directory>] [-n <non_candidate_output_directory>] <file_name>

This script can either run within or without a standard project directory structure:
  -p    If this parameter is provided, this script will run in a "project mode" which needs a standard project directory structure;
  -m    The file name with or without the PATH to a .json files cantains necessary meta information, leave it NULL to find a 'proj_meta.json' file in current WD;

The fillowing parameters is need when it is under a "standalone mode":
  -c    The directory to store the file containing valid reads, leave it NULL to output them in current WD;
  -i    The directory to store the file containing invalid reads for recall, leave it NULL to output them in current WD;

To view the usage information:
  -h    Print usage information and exit.
"""

optlist, args = getopt.getopt(sys.argv[1:], 'phm:c:n:d:')
optdict = dict(optlist)
projWD = os.getcwd()

if ("-h" in sys.argv[1:]) or (len(sys.argv)==1):
    sys.stderr.write(usage)
    sys.exit()

candidate_dir = f"{projWD}/recall/false_split"
discard_dir = f"{projWD}/discard"

if ("-c" in optdict.keys()) and optdict["-c"]:
    candidate_dir = os.path.join(projWD, optdict["-c"])
if ("-n" in optdict.keys()) and optdict["-n"]:
    discard_dir = os.path.join(projWD, optdict["-n"])


# This script will read a sorted list of split reads. 
# Users can either run corresponding commands in advance and run this script under the "standalone" mode,
# or run this script under the "project" mode. 
if "-p" in optdict.keys():  # Project mode
    os.system(f"cat {candidate_dir}/*_true.tsv > {candidate_dir}/split.tsv")
    os.system(f"sort -k 1,1 -T {candidate_dir}/ {candidate_dir}/split.tsv > {candidate_dir}/split.sorted.tsv")
    print(f"[{getDatetime()}] Entries sorted, in directory: '{candidate_dir}'")

    orig_lis = f"{candidate_dir}/split.sorted.tsv"

    candidate_file_name = f"{candidate_dir}/candidate_list.tsv"
    candidate_tooShort_name = f"{candidate_dir}/candidate_tooShort.tsv"
    discarded_file_name = f"{discard_dir}/not_false_split_candidate.tsv"
    onecol_file_name = f"{discard_dir}/one_column.tsv"

else:  # Standalone mode
    orig_lis = args[0]
    print(f"[{getDatetime()}] File: '{orig_lis}' will be processed.")

    candidate_file_name = os.path.join(projWD, candidate_dir, "candidate_list.tsv")
    candidate_tooShort_name = os.path.join(projWD, candidate_dir, "candidate_tooShort.tsv")
    discarded_file_name = os.path.join(projWD, discard_dir, "not_false_split_candidate.tsv")
    onecol_file_name = os.path.join(projWD, discard_dir, "one_column.tsv")



# ================================ Basic Information Loading ====================================
# Loading meta information from a .json file.
if ("-m" in optdict.keys()) and optdict["-m"]:
    meta_json = os.path.join(projWD, optdict["-m"])
else:
    meta_json = f"{projWD}/proj_meta.json"

with open(meta_json) as metaf:
    meta_inf = json.load(metaf)

pattern_basic = {}
for i in meta_inf["UsedAdapter"]:
   pattern_basic[i] = "(?e)(" + meta_inf["AdapterBC"][i] + "){e<=2}"


# ================================= Defining Functions ====================================
# The function to extract essential information from the given entry.
def entryBasicInfExtr(entry):
    entry_lis = entry.strip().split()[0].split("|")
    return (entry_lis[0], int(entry_lis[2]), entry_lis[-1])

# To find a intact BC in an entry.
def adapterAssign4Recall(seq, pdic):
    for bc in pdic.keys():
        bc_hit = regex.search(pdic[bc], seq)
        if bc_hit:
            return (bc, bc_hit.span()[0])

    return None

# Replace the 
def adapterBCassign(entryID, adapterBC):
    ent_lis = entryID.split("|")[:-1]
    ent_lis.append(f"{adapterBC}_falseSplit")
    return "|".join(ent_lis)

def seqLength(seq_entry):
    if len(seq_entry.strip().split())<3:
        return 0
    else: 
        return len(seq_entry.strip().split()[1])

def candidateLenAssess(candidate_entry):
    clen = 0
    for i in candidate_entry:
        clen += seqLength(i)
    return clen>200


# ================================= Defining File Handles ====================================
candidate_file = open(f"{candidate_dir}/candidate_list.tsv", "w")
candidate_tooShort = open(f"{candidate_dir}/candidate_tooShort.tsv", "w")
discarded_file = open(f"{discard_dir}/not_false_split_candidate.tsv", "w")
onecol_file = open(f"{discard_dir}/one_column.tsv", "w")


# ================================= Main loop ====================================
# The main loop of the false split identification process.
# Queues will be established to evaluate which read is potentially false split. 
# A valid queue should has its first element 5' end intact, its last element 3' end intact, other element no intact end.
# All of the elements should come from the same CCS read, and their ID should be continuous.
# When a queue encounters a element with a intact 5' end, this queue should be discarded.
entry_lis = []
case_stat = {
    "Case 1":0, # Case 1：The queue encountered a element with an intact 5' end.
    "Case 2":0, # Case 2：The "first" element of this queue doesn't have an intact 5' end.
    "Case 3":0, # Case 3：This element only has an ID.
    "Case 4":0, # Case 4：The queue encountered an element doesn't have a continuous ID.
    "Case 5":0, # Case 5：Valid queue. 
    "Case 6":0, # Case 6：Valid queue with the sum of its length. 
    "Case 7":0  # Case 7：This element can be added into the current queue.
}

with open(orig_lis) as fil:
    for entry in fil:
        entry_inf = entryBasicInfExtr(entry)

        if entry_inf[-1]=="noBC":
            if len(entry_lis)!=0:
                case_stat["Case 1"] += 1
                for i in entry_lis:
                    discarded_file.write(i)

            entry_lis=[entry]

        else:
            if len(entry_lis)==0:
                case_stat["Case 2"] += 1

                discarded_file.write(entry)
                continue

            entry_sp = entry.split()

            if len(entry_sp)<3:
                case_stat["Case 3"] += 1

                onecol_file.write(entry)
                continue

            # Check if the entry ID is continuous.
            last_inf = entryBasicInfExtr(entry_lis[-1])

            if (entry_inf[0]!=last_inf[0]) or (entry_inf[1]-last_inf[1]!=1):
                case_stat["Case 4"] += 1
                discarded_file.write(entry)

                if len(entry_lis)>0:
                    for i in entry_lis:
                        discarded_file.write(i)
                continue

            resBCassign = adapterAssign4Recall(entry_sp[1], pattern_basic)

            if resBCassign:
                # case_stat["Case 5"] += 1
                entry_lis.append(f"{adapterBCassign(entry.split()[0], resBCassign[0])}\t{entry.split()[1]}\t{entry.split()[-1]}\n")

                if candidateLenAssess(entry_lis):
                    case_stat["Case 5"] += 1
                    for i in entry_lis:
                        candidate_file.write(i)
                else:
                    case_stat["Case 6"] += 1
                    for i in entry_lis:
                        candidate_tooShort.write(i)
                
                entry_lis = []

            else:
                case_stat["Case 7"] += 1
                entry_lis.append(entry)

for i in entry_lis:
    discarded_file.write(i)

candidate_file.close()
candidate_tooShort.close()
discarded_file.close()
onecol_file.close()


# Dump statistic information into a .json file.
print(f"[{getDatetime()}] False split identify process done.")

with open(f"{projWD}/false_split_detect_cases.json", "w") as jf:
    json.dump(case_stat, jf, indent=4)

print(f"[{getDatetime()}] Json file: {projWD}/false_split_detect_cases.json.")