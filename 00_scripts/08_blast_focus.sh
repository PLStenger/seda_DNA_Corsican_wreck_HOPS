#!/bin/bash
#SBATCH --job-name=08_blast
#SBATCH --ntasks=1
#SBATCH -p smp
#SBATCH --mem=250G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL 
#SBATCH --error="/home/plstenge/seda_DNA_Corsican_wreck_HOPS/00_scripts/08_blast.err"
#SBATCH --output="/home/plstenge/seda_DNA_Corsican_wreck_HOPS/00_scripts/08_blast.out"

# Dossier contenant les fichiers fastq.gz
INPUT_DIR="/home/plstenge/seda_DNA_Corsican_wreck_HOPS/06_fastp"
OUTPUT_DIR="/home/plstenge/seda_DNA_Corsican_wreck_HOPS/08_blast"

# Création des dossiers de sortie si besoin
mkdir -p "$OUTPUT_DIR"

 module load conda/4.12.0
 source ~/.bashrc
 conda activate bbduk
 
 REFORMAT=/home/plstenge/bbmap/reformat.sh
 
 cd $INPUT_DIR
 
 gunzip clean_1121_sed8_rep1_dedup_clumpify_fastp_merged.fastq.gz 
 $REFORMAT in=$INPUT_DIR/clean_1121_sed8_rep1_dedup_clumpify_fastp_merged.fastq out=$OUTPUT_DIR/clean_1121_sed8_rep1_dedup_clumpify_fastp_merged.fa

  gunzip clean_1122_sed8_rep2_dedup_clumpify_fastp_merged.fastq.gz 
 $REFORMAT in=$INPUT_DIR/clean_1122_sed8_rep2_dedup_clumpify_fastp_merged.fastq out=$OUTPUT_DIR/clean_1122_sed8_rep2_dedup_clumpify_fastp_merged.fa
 

cd $OUTPUT_DIR

blastn \
  -query clean_1121_sed8_rep1_dedup_clumpify_fastp_merged.fa \
  -db /storage/biodatabanks/ncbi/NT/current/flat/nt \
  -out clean_1121_sed8_rep1_dedup_clumpify_fastp_merged.fa.tsv \
  -outfmt 6 \
  -max_target_seqs 10 \
  -num_threads 4

blastn \
  -query clean_1122_sed8_rep2_dedup_clumpify_fastp_merged.fa \
  -db /storage/biodatabanks/ncbi/NT/current/flat/nt \
  -out clean_1122_sed8_rep2_dedup_clumpify_fastp_merged.fa.tsv \
  -outfmt 6 \
  -max_target_seqs 10 \
  -num_threads 4
