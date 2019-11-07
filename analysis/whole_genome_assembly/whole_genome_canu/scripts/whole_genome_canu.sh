#!/bin/bash

CANU=/well/longread/users/akl399/bin/canu-1.9/Linux-amd64/bin/canu

FASTQ=/well/longread/projects/nanopore/data/WTON000393/fastq-trimmed-all.fastq.gz

OUTPUT_FOLDER=data/ONT

module load gcc/5.4.0
module load java/1.8.0_latest

$CANU -p canu-ONT-WGA -d $OUTPUT_FOLDER genomeSize=3.2g 'gridOptions=-q himem.qh -P todd.prjc -S /bin/bash -N canu-ONT-WGA -e canu_ONT_WGA_error.log' 'gridEngineResourceOption=-pe shmem THREADS' -nanopore-raw $FASTQ
