#!/bin/bash
#SBATCH --job-name=11_Bowtie2_indexation
#SBATCH --ntasks=1
#SBATCH -p smp
#SBATCH --mem=1000G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL 
#SBATCH --error="/home/plstenge/seda_DNA_Corsican_wreck_HOPS/00_scripts/11_Bowtie2_indexation.err"
#SBATCH --output="/home/plstenge/seda_DNA_Corsican_wreck_HOPS/00_scripts/11_Bowtie2_indexation.out"

# Activation environnement
module load conda/4.12.0
source ~/.bashrc
conda activate metaDMG

bowtie2-build /storage/biodatabanks/ncbi/NT/current/fasta/All/all.fasta /home/plstenge/refmetagenome
