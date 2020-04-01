#!/bin/bash
#SBATCH
#SBATCH --job-name=Phage
#SBATCH --time=72:00:00
#SBATCH --partition=lrgmem
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mail-type=end
#SBATCH --mail-user=dmonaco1@jhu.edu
#SBATCH --output=/data/hlarman1/PhIPdb/Studies/Daniel/Phageome/Blast/log/%J-%x.out


#
#---------------------------------------------------------------------
# SLURM job script to run serial R
#---------------------------------------------------------------------

ml R/3.6.1

ml

export PATH=/data/hlarman1/PhIPdb/Software/Alignment/BLAST+/ncbi-blast-2.8.1+/bin:"$PATH"
export PATH=/data/hlarman1/PhIPdb/Software/texlive/2018/bin/x86_64-linux:"$PATH"
export PATH=/data/hlarman1/PhIPdb/Software/BLAST+/ncbi-blast-2.8.1+/bin:"$PATH"

Rscript /data/hlarman1/PhIPdb/Studies/Daniel/Phageome/Blast/phage_blast_test.R
