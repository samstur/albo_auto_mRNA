#!/bin/bash
#SBATCH --job-name=htseq_count --output=%x.%j.out
#SBATCH --mail-type=END,FAIL --mail-user=sls366@georgetown.edu
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --time=24:00:00
#SBATCH --mem=4G

#-----------------------------------------------------------------------------#
# This script uses htseq to count and assign reads to genes #
#-----------------------------------------------------------------------------#

source activate conda_env

#- Set variables ----------------------------------------------------------------#

bam_dir=/home/sls366/albo_biting/mRNA/bam_dir
count_dir=/home/sls366/albo_biting/mRNA/counts_dir
htseq=/home/sls366/.conda/envs/conda_env/bin/htseq-count
ref=/home/sls366/albo_biting/genome_files/aedes_albopictus_AalbF3.gff3

#- RUN fastqc ----------------------------------------------------------------#

files=(${bam_dir}/*_Aligned.out.bam)
for file in ${files[@]}
do
base=`basename ${file} _Aligned.out.bam`
${htseq} -f bam -r pos -s yes -t exon -i gene ${bam_dir}/${base}_Aligned.out.bam ${ref} > ${count_dir}/${base}_htseqCount

done

## -s yes indicates that our sequencing is forward stranded

#- FIN -----------------------------------------------------------------------#
