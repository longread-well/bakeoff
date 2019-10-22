#!/bin/bash

# Specify a job name
#$ -N Whole_gemome_canu

# Project name and target queue
#$ -P todd.prjc
#$ -q himem.qh

# Run the job in the current working directory
#$ -cwd -j y

# Log locations which are relative to the current
# working directory of the submission
#$ -o qsub_output.log
#$ -e qsub_error.log

# Parallel environemnt settings
#  For more information on these please see the wiki
#  Allowed settings:
#   shmem
#   mpi
#   node_mpi
#   ramdisk
#$ -pe shmem 6

# Some useful data about the job to help with debugging
echo "------------------------------------------------"
echo "SGE Job ID: $JOB_ID"
echo "SGE Task ID: $SGE_TASK_ID"
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
echo "------------------------------------------------"

# Begin writing your script here

FASTQ=/well/longread/projects/nanopore/data/WTON000393/fastq-trimmed-all.fastq.gz
OUTPUT_FOLDER=data/Canu

module load gcc/5.4.0
module load java/1.8.0_latest

/well/ont/apps/canu-1.8/Linux-amd64/bin/canu -d $OUTPUT_FOLDER -p canu genomeSize=3.2g 'useGrid=false' -nanopore-raw $FASTQ
echo "Assembly complete"

# End of job script
