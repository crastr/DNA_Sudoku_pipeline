
Welcome to the **dePoP** wiki!

## Overview
**dePoP** is a pipeline designed for analysis of NGS-sequences of pooled samples aimed on identification of rare Nucleotide Variants (NV) carriers in cases of usage of an overlapping pooling scheme.

**dePoP** automatizes de-multiplexing, trimming, mapping, snp-calling of raw NGS-reads and performs de-pooling, i.e. identification of carrying NV individual samples. De-pooling process is performed by ***s-dePooler*** - a novel java application.

**dePoP** is adapted for work in Linux OS.
## Requirements
**dePoP** is designed for work in Linux OS. Perl interpreter is needed. 
The pipeline itself requires no installation. But the following programs must be registered in the environment:
* [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) ver 2.2 or newer 
* [Cutadapt](http://cutadapt.readthedocs.io/en/stable/guide.html)  ver 1.8 or newer
* [samtools](https://sourceforge.net/projects/samtools/files/samtools/) and [bcftools](https://sourceforge.net/projects/samtools/files/samtools/) ver 1.3 or newer.

The following tools are included in repository and alternatively can be installed separately and included in pipeline by defined options:
* [Genome Analysis Toolkit (GATK)](https://software.broadinstitute.org/gatk/) .
* ***s-dePooler***
