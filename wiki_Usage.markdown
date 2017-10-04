### Files needed
1) Reads in solexa/454 format
2) Reference in fasta format
3) Bed-file containing regions of interest for reference
# ??? ïðèìåð ôàéëà(ìîæíî ññûëêîé)
4) A list of files containing comma-separated sequences to trim from 5’end.
5) A list of files containing comma-separated sequences to trim from 3’end.
6) Fasta-file containing demultiplication primers. 
7) A file containing adapters to cut in fasta format 

### Operation
1) Create a BowTie2 index of the reference file using the following command
```
bowtie2-build -f [Your reference fasta file] [Your desired output prefix]
```
2) Create a pooling table of the following architecture:
# ?????

3) Create a cfg.file containing all  information necessary for the run.
4) Run the pipeline like that: `perl *.pl --cfg [path to cfg file]`

###### Path to cfg file
|Parameter|Type|Description|
|:---|:---:|:---|
|--cfg | string | Path to cfg-file. Can be overridden with command line parameters.|

###### Procedures
|Parameter|Type|Description|
|:---|:---:|:---|
|--dmpf|[string]|Path to folder containing demultiplexed reads. If set - script starts with “Reads mapping” step.|
|--bam|[string]|Path to merged bam-file, reads of different pools must be marked by ?? tag. If set - script starts with “SNP-calling” step.|
|--vcf|[string]|Path to VCF file. If set - script starts with depooling step.|
###### Run_parameters
|Parameter|Type|Description|
|:---|:---:|:---|
|--prefix|[string]|Prefix added to all generated files, “run” by default.|
|--wf|[string]|Work folder.|
|--gm|[flag]|Graphic monitor on depooling process.|
|--output|[string]|Output directory|
###### Data
|Parameter|Type|Description|
|:---|:---:|:---|
|--fq|[string]|Path to the fastq read file.|
|--fa|[string]|Path to fasta containing 454 reads|
|--q|[string]|Path to quality file for 454 reads|
|--scheme|[string]|Path to file containig pooling scheme.|
|--demult|[string]|Path to fasta-file containing demultiplication primers. The beginning of primer sequence header will be used for pool naming|
|--trseq5|[string]|A list of files containing comma-separated sequences to trim from 5’end.|
|--trseq3|[string]|A list of files containing comma-separated sequences to trim from 3’end.|
|--trer|[float]|Value of trimming error tolerance. Value ] 0  Default 0.2.|
|--trqual|[integer]|Value of trimming quality  (Cutadapt -q) Default 25.|
|--readlen|[integer]|Read length filter (cutadapt -m) Default 40.|
|--bt2index|[string]|Path to the the bowtie index file|
|--reference|[string]|Path to the reference fasta file.|
|--regions|[string]|Path to bed file, containing coordinates of regions of interest*|
|--org_ploidy|[integer]|Ploidy of a single sample (ploidy of an organism)|
|--pool_size|[integer]|Number of samples in a pool.|
|--T|[string]|Options for SNV-calling. Possible values: HaplotypeCaller, or ???|
|--gt_mode [string]|[string]|Options for SNV-calling gt_mode. . Possible values: “???” or ???|

###### Miscelanious
|Parameter|Type|Description|
|:---|:---:|:---|
|--n| [int] |max number of threads (2 by default)|
|--gatk| [string] |Custom GATK location. If not set the included version will be used.|
|--454tofq| [string] |Path to a 454-to-fastq converting script|