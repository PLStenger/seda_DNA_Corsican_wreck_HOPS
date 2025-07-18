#!/bin/bash
#SBATCH --job-name=08_krona
#SBATCH --ntasks=1
#SBATCH -p smp
#SBATCH --mem=700G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL 
#SBATCH --error="/home/plstenge/seda_DNA_Corsican_wreck_HOPS/00_scripts/08_krona.err"
#SBATCH --output="/home/plstenge/seda_DNA_Corsican_wreck_HOPS/00_scripts/08_krona.out"

module load conda/4.12.0
source ~/.bashrc
conda activate krona

# For installation
# https://telatin.github.io/microbiome-bioinformatics/Kraken-to-Krona/ :
# conda install -c bioconda krona
# then run ktUpdateTaxonomy.sh

# Répertoires
IN_DIR="/home/plstenge/seda_DNA_Corsican_wreck_HOPS/07_kraken2"
OUT_DIR="/home/plstenge/seda_DNA_Corsican_wreck_HOPS/08_krona"

mkdir -p "$OUT_DIR"

cd $IN_DIR

ktImportTaxonomy -t 5 -m 3 -o $OUT_DIR/krona.html $IN_DIR/clean_1129_sed6_rep1_dedup_clumpify_merged.report
echo "Terminé : clean_1129_sed6_rep1_dedup_clumpify_merged.report"

ktImportTaxonomy -t 5 -m 3 -o $OUT_DIR/multi-krona.html $IN_DIR/*.report 
echo "Terminé : all"
