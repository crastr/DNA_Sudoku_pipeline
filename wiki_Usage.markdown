## Launch
perl dna_sudoku_devel.pl [options]
## Process
Pipeline **dePoP** comprises several sequential stages:
 - raw reads demultiplexing and trimming
 - reads mapping
 - snp-calling
 - depooling.
Analysis can be started from any of these stages that is defined by a provided files.

## Files needed
(Bold parameters are mandatory for a corresponding stage and following stages. )

(Parameters marked with * are key for a corresponding stage. If such parameters is defined the process will start form this stage. )
### Pooling scheme 
|Parameter| Description |
|:--|:--|
| **-scheme** | File containing pooling scheme in a tab-separated format. (See below.) |
### Stages
#### Raw reads demultiplexing and trimming
|Parameter| Description |
|:--|:--|
| **-fq\*** | Fastq-file with raw reads. |
| **-demult** | Fasta-file containing demultiplication sequences. Each demultiplexing sequence corresponds to a pool. The first token before any white-space in the header of the sequence must be the pool identificator in the Pooling scheme file. |
| -trseq5 | Fasta-file containing adapter sequences to trim from 5’end. If sequences require sequential adapters trimming, this parameters can be specified more than once or different files can be separated by a comma. |
| -trseq3 | Fasta-file containing adapter sequences to trim from 3’end. If sequences require sequential adapters trimming, this parameters can be specified more than once or different files can be separated by a comma. |

#### Reads mapping
|Parameter| Description |
|:--|:--|
| -dmpf\* | Path to folder containing demultiplexed reads. Fastq-files should be named as pool identificators in Pool  Scheme file. File expansion '*.fq'|
| **-bt2index** | Path to the the bowtie index file |

#### Snp-calling
|Parameter| Description |
|:--|:--|
| -bam\* | Path to a single bam-file containing all mapped reads. Reads of different pools must be marked by 'RG'- tag corresponding to the pool identificator as defined in Pool  Scheme file. |
| **-reference** | Path to the reference fasta file.  fasta-file should be indexed.|
| **- regions** | Path to bed file, containing coordinates of regions of interest. |
#### Depooling
|Parameter| Description |
|:--|:--|
| -vcf\* | Path to VCF-file. File should contains all Pool identificators as Samples. File must be indexed. |
| - regions | Path to bed file, containing coordinates of regions of interest. |

### Configuration file
Configuration file is a simple text file. Each line comprises an option and its value separated by a white-space. Values defined in a configuration file will be overridden by parameters of command line.
|Parameter| Description |
|--|--|
| -cfg | Path to configuration file. |

### dePoP options
|Parameter|Type|Description|
|:---|:---:|:---|
|-n| [integer] |max number of threads (2 by default)|
|-prefix|[string]|Prefix added to all generated files, “run” by default.|
|-output|[string]|Output directory|
|-org_ploidy|[integer]|Ploidy of a single sample (ploidy of an organism) (2 by default)|
| -gatk	| [string]	| Custom GATK location. If not set the included version will be used. |
### Process options
### Trimming/reads filtering options
|Parameter|Type|Description|
|:---|:---:|:---|
|-trer|[float]|Value of trimming error tolerance. Value ] 0  Default 0.2.|
|-trqual|[integer]|Value of trimming quality  (Cutadapt -q) Default 25.|
|-readlen|[integer]|Read length filter (cutadapt -m) Default 40.|

### Depooling options
|Parameter|Type|Description|
|:---|:---:|:---|
| -seqerr | [float] | Sequencing error rate. Value (0; 0.2]|
| -conf | [float] | Confidence interval inside Allele-Pool distribution. Value (0.6; 0.95] |
| -mp | [float] | Mixing accuracy Value [1; 10]|
| -time_limit | [integer] | Time in millisecunds for combination check. |
| -cmb_limit | [integer] |  Power of number combination that will be checked before interapted. |
| -tries | [integer] | Maximal number of variants distributions of alleles between pools. |
| -success | [integer] | Number of successfully depooled tries. |
| -noIndels | [flag] | No insertion or deletions polymorphisms. |
|-gm|[flag]|Graphic monitor on depooling process.|

