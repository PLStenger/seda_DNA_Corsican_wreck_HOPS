#!/bin/bash
#SBATCH --job-name=06_kraken_db
#SBATCH --ntasks=1
#SBATCH -p smp
#SBATCH --mem=250G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL 
#SBATCH --error="/home/plstenge/seda_DNA_Corsican_wreck_HOPS/00_scripts/06_kraken_db.err"
#SBATCH --output="/home/plstenge/seda_DNA_Corsican_wreck_HOPS/00_scripts/06_kraken_db.out"

module load conda/4.12.0
source ~/.bashrc
conda activate kraken2

cd /home/plstenge/nt_kraken2_db

kraken2-build --add-to-library /storage/biodatabanks/ncbi/NT/ncbi_blast_nt_2024-8-24/fasta/All/all.fasta --db /home/plstenge/nt_kraken2_db

echo "add-to-library ok"

kraken2-build --download-taxonomy --db /home/plstenge/nt_kraken2_db

echo "download-taxonomy ok"

kraken2-build --build --db /home/plstenge/nt_kraken2_db --threads 36

echo "base de donnée construite"
