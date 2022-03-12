#!/bin/bash
#SBATCH --job-name=fastqc --output=%x.%j.out
#SBATCH --mail-type=END,FAIL --mail-user=sls366@georgetown.edu
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --time=12:00:00
#SBATCH --mem=1G

#-----------------------------------------------------------------------------#
# This script gives quality control of fastq files #
#-----------------------------------------------------------------------------#


#- RUN fastqc ----------------------------------------------------------------#

/home/sls366/FastQC/fastqc -o /home/sls366/albo_biting/mRNA/fastqc_dir /home/sls366/albo_biting/mRNA/trim_dir/*PE.fastq.gz

#- FIN -----------------------------------------------------------------------#
