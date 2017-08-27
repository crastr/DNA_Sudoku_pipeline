#!/bin/sh
prefix=$1
sam=$prefix'.sam'
report=$prefix'_report.txt'
uns_bam=$prefix'_unsorted.bam'
bam=$prefix'.bam'

if [ -f $bam ]; then rm -f $bam; fi
if [ -f $bam.bai ]; then rm -f $bam.bai; fi
if [ -f $report ]; then rm -f $report; fi

bowtie2-align-s -x $2 -U $3 --rg-id SM --rg 'SM:'$4 --threads $5 -S $sam 2>> $report 

samtools view -bS $sam > $uns_bam
samtools 'sort' $uns_bam -o $bam
samtools index $bam

rm $sam $uns_bam
#chmod -w $report $bam $bam.bai
