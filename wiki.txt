# Manual
## Review
A pipeline for analysis of binned samples of different origin. Complete pipeline transforms reads in 454/Solexa formats into tables, showing correspondence between polymorphisms and individual samples(patients). 
## Installation
DNA_Sudoku needs the followng programs in your $PATH
* Bowtie2 ver X or newer http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
* Cutadapt http://cutadapt.readthedocs.io/en/stable/guide.html
* Samtools http://samtools.sourceforge.net/
* GATK.jar (optional) ???
## Usage
### Files needed
* Reads in solexa/454 format
* Reference in fasta format
* Bed-file containing regions of interest for reference
* A list of files containing comma-separated sequences to trim from 5’end.
* A list of files containing comma-separated sequences to trim from 3’end.
* Fasta-file containing demultiplication primers. 
* File containing all the primers ???
* File containing all the primers ???
### Operation
1) Create a BowTie2 index of the reference file using the following command
|| ???code here!!! ||
2) Create a pooling table of the following architecture
3) Create a cfg.file containing all necessary information for the run.
4) Run the pipelinr like: perl *.pl --cfg <path to cfg file>

--cfg|<string>|Path to cfg-file. Can be overridden with command line parameters.|
-----configuration parameters-----|
|--n|<int>|max number of threads (2 by default)|
|--gatk|<string>|Custom GATK location. If not set the included version will be used.|
|--454tofq|<string>|Path to a 454-to-fastq converting script|
----Procedures----
|--fq|<string>|Path to the fastq read file.|
|--dmpf|<string>|Path to folder containing demultiplexed reads. If set - script starts with “Reads mapping” step.|
|--bam|<string>|Path to merged bam-file, reads of different pools must be marked by ?? tag. If set - script starts with “SNP-calling” step.|
|--vcf|<string>|Path to VCF file. If set - script starts with depooling step.|
# Run_parameters----
|--prefix|<string>|Prefix added to all generated files, “run” by default.|
|--wf|<string>|Work folder.|
|--gm|<boolean>|Graphic monitor on depooling process.|
|--output|<string>|Output directory|
---Data---
|--fa|<string>|Path to fasta containing 454 reads|
|--q|<string>|Path to quality file for 454 reads|
|--scheme|<string>|Path to file containig pooling scheme.|
|--demult|<string>|Path to fasta-file containing demultiplication primers. Pool identificator should be the first token in a sequence header.|
|--trseq5|<string>|A list of files containing comma-separated sequences to trim from 5’end.|
|--trseq3|<string>|A list of files containing comma-separated sequences to trim from 3’end.|
|--trer|<float>|Value of trimming error tolerance. Value > 0  Default 0.2.|
|--trqual|<integer>|Value of trimming quality  (Cutadapt -q) Default 25.|
|--readlen|<integer>|Read length filter (cutadapt -m) Default 40.|
|--bt2index|<string>|Path to the the bowtie index file|
|--reference|<string>|Path to the reference fasta file.|
|--regions|<string>|Path to bed file, containing coordinates of regions of interest*|
|--org_ploidy|<integer>|Ploidy of a single sample (ploidy of an organism)|
|--pool_size|<integer>|Number of samples in a pool.|
|--T|<string>|Options for SNV-calling. Possible values: HaplotypeCaller, or ???|
|--gt_mode <string>|<string>|Options for SNVyou-calling gt_mode. . Possible values: “???” or ???|
|--454tofq|<STRING>|Path to a custom 454-to-fastq converting script|

Bin:|*.pl - скрипт|SudokuSolver.jar|GATK.jar (optional)|
