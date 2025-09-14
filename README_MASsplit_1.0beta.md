# Workflow for Splitting MAS-PAIso-Seq(2) Data

This workflow contains several scripts to extract, split and recall transcript reads from an unaligned `.bam` file containing converted MAS-PAIso-seq(2) CCS reads.
This is a beta version of the MAS-PAIso-seq(2) data splitting workflow.

If more funcyion is needed while using this working flow, please contact:

``` text
jiazheng@genetics.ac.cn
```

## Software Dependencies

This splitting workflow depends on these softwares and Python third-party libraries:

``` text
# Software & Script
fastq-splitter.pl (http://kirill-kryukov.com/study/tools/fastq-splitter)
ParaFly (http://parafly.sourceforge.net/)

# Python third-party library
regex  (2024.9.11)
pysam  (0.22.1)
```

## Preparing the meta-information file

A meta-information file is needed when running the following scripts in the current splitting workflow:

- `split_MASseq_v1.0a.py`
- `recall_MASseq_v1.0a.py`
- `false_split_detect_v1.0a.py`
- `convert_tsv2fqgz_v1.0a.py`

A general template of the `.json` file is:

``` json
{
    "Project": "<Project_name>",
    "Approach": "MAS-PAIso-seq(2)",
    "Description": "<Project_description>",  

    "AdapterBC": {
        "BC1": "TGCTATCTGAGATACT",
        "BC2": "GAGTCTCGATATACTA",
        "BC3": "ACGTGCTCTATAGAGA",
        "BC4": "TATCAGCACGACATGC",
        "BC5": "GCTCTCACGATATCAG",
        "BC6": "TATATGCTCTGTGTGA",
        "BC7": "ATGATGTGCTACATCT",
        "BC11": "GCGCACGCACTACAGA",
        "BC12": "ACACTGACGTCGCGAC",
        "BC13": "CGTCTATATACGTATA",
        "BC14": "ATAGAGACTCAGAGCT",
        "BC15": "TAGATGCGAGAGTAGA"
    },
    "Adapter3GeneralSeq": "GTACTCTGCGTTGATACCACTGCTT",

    "UsedAdapter": ["BC<x>", "BC<y>"],
    
    "Adapter2Sample": {
        "BC<x>": "<sample_ID_x>", 
        "BC<y>": "<sample_ID_y>"
    }
}
```

This meta-information file contains the following essential keys:

- `AdapterBC <dic>`: The standard sequence of 3' adapter barcode, these barcodes will be used to identify which sample the read comes from.
- `Adapter3GeneralSeq <str>`: The sequence following the barcode of a 3'adapter, which is common among all the 3' adapters.
- These essential keys are necessary for:
  - `split_MASseq_v1.0b.py`
  - `recall_MASseq_v1.0b.py`
  - `false_split_detect_v1.0b.py`

and a flexible key:

- `UsedAdapter <list>`: 3' adapters used in the sequencing project, which is denoted by barcode IDs consistent with the keys of `AdapterBC`.
- This key is necessary for:
  - `split_MASseq_v1.0b.py`
  - `recall_MASseq_v1.0b.py`

After all, there are some optional keys to improve the using experience of this workflow:

- `Adapter2Sample <dic>`: The corresponding relationship of barcode IDs (as keys) and sample names (as values).
- When this key is provided, output files of this script will named after their samples instead of their barcodes.

Other keys are provided soely to improve the readability of this file.

This `.json` file can have any name when the workflow run under the `standalone` mode. When it is used under a `project` mode, its name has to be `proj_meta.json`.

## Running the workflow

This workflow can be run either under a `project mode` or a `standalone mode`.
When this workflow run under a `standalone mode`, it means scripts will be separatly run without a standard directory structure.
When this workflow run under a `project mode`, a proper dierctory should be established, with scripts and data files correctly placed.
A `.json` file containing necessary meta information should always prepared to run this workflow.

### Option 1. Running under `standalone mode`

When running this workflow under a `standalone mode`, scripts will run separately, thus the workflow can be run step-by-step manually, with a directory structure more flexible than `project` mode.

#### Step 1.1. Read Extracting

Note that current PacBio sequencing result is usually provided in a converted `unaligned BAM` file rather than a raw `PacBio BAM` format, thus CCS reads can be directly extracted. If your data is in a raw `PacBio BAM` format, please refer to our previous workflow to conduct format converting.

Script `extr_MASseq_v1.0b.py` is used to extract CCS reads from the `.bam` file, which will identify the `.bam` file in the `hifi_reads/` directory and extract the name, sequence, quality and pass number of a CSS read, output these information to a `.fastq` file:

``` bash
python extr_MASseq_v1.0b.py -o <output_directory> -f <output_filename> [<PATH>/]<file_name>.bam
```

For more help information, please run `python extr_MASseq_v1.0b.py -h`.

#### Step 1.2. CCS Read Splitting

Our current MAS-PAIso-seq(2) will product libraries have following structure:
![[MAS-PAIso-seq(2)_lib_struc.png]]

During the library construction process, a 3' adapter sequence and a 5' UMI is firstly added to the 3' end of the transcript (with their common sequence noted as `SigRc` and `SigF`, respectively), then a 5' MAS primer and a 3' MAS primer was added.

The CCS reads will be split upon the `SigRc` sequence. It is identified as a "valid" read if:

- It has an intact `SigF` sequence in its first 50 nt, which means it has an mostly intact 5' end;
- It has an identifiable `Adapter Barcode` in its last 100 nt, which means it has an intact 3' end;
- A `UMI` pattern can be find right after the `SigF` sequence.

Script `split_MASseq_v1.0b.py` is used to split CCS reads and validate split results:

``` bash
python split_MASseq_v1.0b.py -m <meta_information_json> -v <valid_output_directory> -i <invalid_output_directory> [<PATH>/]<file_name>.fastq
```

The output format of this script is a `.tsv` table containing following three columns, which can be conveniently converted to `.fastq` format:

``` text
<read_ID>  <read_sequence>  <read_quality_string>
```

The `<read_ID>` column has serveral fields separated by "`|`":

``` text
@<CCS_read_ID>|<pass_num>|<seq_num>|<barcode_ID>|<UMI_seq>
```

- `CCS_read_ID`: the ID of the original CCS read where this split read came from.
- `pass_num`: the pass number of the original CCS read.
- `seq_num`: if this split read is the `n` th read splitted from its original CCS read, this value is set as n.
- `barcode_ID`: the ID of the barcode detected from the split read.
- `UMI_seq`: the UMI sequence in the 5' end of the split read.

For more help information, please run `python extr_MASseq_v1.0b.py -h`.

#### Step 1.3. Read Recalling

Some potentially valid reads may detected as "invalid" in the upstream CCS read splitting step. These reads can be recalled using a slightly altered validating protocol, which was implemented in script `recall_MASseq_v1.0b.py`:

``` bash
python recall_MASseq_v1.0b.py -m <meta_information_json> -r <recall_files_directory> -v <valid_output_directory> -d <discarded_output_directory> <file_name>
```

If this script run under the `standalone` mode, its input files should be prepared in advance:

``` bash
cat <invalid_output_directory>/*.err.tsv > <recall_files_directory>/<file_name>.err.tsv
cat <invalid_output_directory>/*.deg.tsv > <recall_files_directory>/<file_name>.deg.tsv
cat <invalid_output_directory>/*.noBC.tsv > <recall_files_directory>/<file_name>.noBC.tsv
```

For more help information, please run `python recall_MASseq_v1.0b.py -h`.

#### Step 1.4. (Optional) False Split Detecting

There is a potential case that some transcripts may contain the `SigRc` sequence in some species, which lead to a “false split” event that split a valid read into two invalid reads. Thus the detection of false splits will consider invalid reads that cannot be recalled.
False split detect protocol is implenemted in script `false_split_detect_v1.0b.py`:

``` bash
python split_MASseq_<version>.py -m <meta_information_file> -c <candidate_output_directory> -n <non_candidate_output_directory> <file_name_sorted>
```

This script will create a file names `candidate_list.tsv` containing false split candidates. There are few candidates detected in our previous works. If a large amount of candidates are detected in your samples, please contact us for a script to recall them.

If this script run under the `standalone` mode, its input file should be prepared in advance:

``` bash
cat *_true.tsv > <PATH>/<file_name>  # "*_true.tsv" recfers to the files containing reads that cannot recalled by the pervious script "".
sort -k 1,1 <PATH>/<file_name> > <PATH>/<file_name_sorted>
```

For more help information, please run `python recall_MASseq_v1.0b.py -h`.

#### 1.5 Splitting Result Reformatting

If the workflow is run under a `standalone` mode, following reformatting commands are recommended to extract and reformat the "valid" reads information of each sample from the above mentioned `.tsv` files to `.fastq.gz` files for further analysis:

``` bash
awk -F '|' '$4=="BC<num>"' <valid_output_directory>/*.tsv | awk '{print "@"$1"\n"$2"\n""+""\n"$3}'| gzip -nc > <sample_name>.fastq.gz &
```

### Option2. Running under `project mode`

A more integrated mode called `project mode`, is more recommended to run this splitting workflow. A few steps are needed to establish a workflow under `project` mode.

A workflow under `project` mode is typically established in the `PacBio` sequencing data directory. Such directories usually have its structure as below:

``` text
.
└── cell-<x>/
    └── hifi_reads/
        ├── <filename>.bam
        ├── <filename>.bam.pbi
        └── <filename>.bam.xml
```

The `bash` script `run_project_mode_v1.0b.sh` and the meta-information file `proj_meta.json` should be copied to the `cell-<x>` directory, while a directory `scripts/` containing these scripts should be created manually:

- `extr_MASseq_v1.0b.py` - provided in this repository.
- `split_MASseq_v1.0b.py` - provided in this repository.
- `recall_MASseq_v1.0b.py` - provided in this repository.
- `false_split_detect_v1.0b.py` - provided in this repository.
- `convert_tsv2fqgz_v1.0b.py` - provided in this repository.
- `fastq-splitter.pl` - can be download from <http://kirill-kryukov.com/study/tools/fastq-splitter>

Then create or modify the `proj_meta.json` according to the experiment design, the "flexibile" key `UsedAdapter` should be a `list` object containing the 3' adapter barcodes' ID used in the experiment, and the "optional" key `Adapter2Sample` are suggested to be a `dictionary` object with its keys as the 3' adapter barcodes' ID and its values as sample names.

After the above mentioned perparation process, the `cell-<x>/` will have a structure as:

``` text
cell-<x>/
├── hifi_reads/
│   ├── <filename>.bam
│   ├── <filename>.bam.pbi
│   └── <filename>.bam.xml
├── run_project_mode_v1.0b.sh
├── proj_meta.json
└── scripts/
    ├── extr_MASseq_v1.0b.py
    ├── fastq-splitter.pl
    ├── split_MASseq_v1.0b.py
    ├── recall_MASseq_v1.0b.py
    ├── false_split_detect_v1.0b.py
    └── convert_tsv2fqgz_v1.0b.py
```

Finally, run this command under `cell-<x>/` to run the splitting workflow under the `project` mode:

``` bash
bash run_project_mode_v1.0b.sh
```

This workflow will automatically establish a standard "project" directory structure and conduct the splitting and recalling process. The resuliting `.fastq.gz` files can be found in the `split_result/` directory.
