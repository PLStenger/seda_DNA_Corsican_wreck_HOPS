#!/bin/bash
#SBATCH --job-name=09_HOPS_homemade
#SBATCH --ntasks=1
#SBATCH -p smp
#SBATCH --mem=1000G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL 
#SBATCH --error="/home/plstenge/seda_DNA_Corsican_wreck_HOPS/00_scripts/09_HOPS_homemade.err"
#SBATCH --output="/home/plstenge/seda_DNA_Corsican_wreck_HOPS/00_scripts/09_HOPS_homemade.out"

module load conda/4.12.0
source ~/.bashrc
conda activate hops

# Dossier contenant les fichiers Kraken2
KRAKEN_DIR="/home/plstenge/seda_DNA_Corsican_wreck_HOPS/07_kraken2"
# Dossier contenant les FASTQ
FASTQ_DIR="/home/plstenge/seda_DNA_Corsican_wreck_HOPS/06_fastp"
# Répertoire de sortie pour les profils de dommage
DAMAGE_DIR="/home/plstenge/seda_DNA_Corsican_wreck_HOPS/08_damage"
mkdir -p $DAMAGE_DIR

LOGFILE="/home/plstenge/seda_DNA_Corsican_wreck_HOPS/00_scripts/09_HOPS_homemade_$(date +%Y%m%d_%H%M%S).txt"
echo "Logging to $LOGFILE"
echo "Script started at $(date)" | tee -a "$LOGFILE"

#declare -A TAXONS=(
#    ["Mammals"]="40674"
#    ["Homo sapiens"]="9606"
#    ["Magnoliopsida"]="3398"
#    ["Vitaceae"]="3602"
#    ["Vitis vinifera"]="29760"
#    ["Bacteria"]="2"
#    ["Actinopterygii"]="7898"
#)

# Liste de groupes taxonomiques à analyser
# Tableau associatif : taxon (clef) -> tax_id (valeur) -> génome de référence (chemin)
declare -A TAXONS=(
    ["Homo_sapiens"]="9606:/storage/biodatabanks/ucsc/genomes/hg19/Homo_sapiens-hg19_2012-9-19/fasta/all.fasta"
    ["Vitis_vinifera"]="29760:/storage/groups/gdec/shared_paleo/genomes_REF/12Xv2_grapevine_genome_assembly.fa"
)


touch "$LOGFILE"
echo "Logging to $LOGFILE" | tee -a "$LOGFILE"

for KRAKEN_FILE in "$KRAKEN_DIR"/*.kraken; do
    KRAKEN_BASE=$(basename "$KRAKEN_FILE" .kraken)
    echo -e "\n==== Processing file: $KRAKEN_FILE ====" | tee -a "$LOGFILE"

    # Extraire le préfixe FASTQ avant _dedup_clumpify_unmerged ou _dedup_clumpify_merged
    PREFIX=$(echo "$KRAKEN_BASE" | sed -E 's/_dedup_clumpify_(un)?merged$//')
    echo "Prefix for FASTQ search: $PREFIX" | tee -a "$LOGFILE"

    R1_FILE="${FASTQ_DIR}/${PREFIX}_dedup_clumpify_fastp_R1.fastq.gz"
    R2_FILE="${FASTQ_DIR}/${PREFIX}_dedup_clumpify_fastp_R2.fastq.gz"
    MERGED_FILE="${FASTQ_DIR}/${PREFIX}_dedup_clumpify_fastp_merged.fastq.gz"

    echo "Looking for R1: $R1_FILE" | tee -a "$LOGFILE"
    echo "Looking for R2: $R2_FILE" | tee -a "$LOGFILE"
    echo "Looking for merged: $MERGED_FILE" | tee -a "$LOGFILE"

    for GROUP in "${!TAXONS[@]}"; do
        TAX_ID="${TAXONS[$GROUP]%:*}"
        REF_FASTA="${TAXONS[$GROUP]#*:}"
        echo "--- Processing taxon: $GROUP (taxid: $TAX_ID) ---" | tee -a "$LOGFILE"

        OUT_R1="${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_R1.fastq"
        OUT_R2="${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_R2.fastq"
        OUT_MERGED="${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.fastq"

        # Extraction pour reads appariés (unmerged)
        if [[ -f "$R1_FILE" && -f "$R2_FILE" ]]; then
            echo "Extracting reads for $GROUP (unmerged pairs)..." | tee -a "$LOGFILE"
            python3 /home/plstenge/KrakenTools/extract_kraken_reads.py \
                -k "$KRAKEN_FILE" -s "$R1_FILE" -s2 "$R2_FILE" -t "$TAX_ID" \
                -o "$OUT_R1" -o2 "$OUT_R2" --fastq-output 2>>"$LOGFILE"

            if [[ -f "$OUT_R1" && -f "$OUT_R2" ]]; then
                echo "Aligning with bwa aln/sampe (unmerged pairs)..." | tee -a "$LOGFILE"
                bwa aln -t 4 "$REF_FASTA" "$OUT_R1" > "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_R1.sai" 2>>"$LOGFILE"
                bwa aln -t 4 "$REF_FASTA" "$OUT_R2" > "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_R2.sai" 2>>"$LOGFILE"
                bwa sampe "$REF_FASTA" "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_R1.sai" "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_R2.sai" "$OUT_R1" "$OUT_R2" > "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}.sam" 2>>"$LOGFILE"
                samtools view -bS "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}.sam" > "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}.bam" 2>>"$LOGFILE"
                samtools sort -o "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}.sorted.bam" "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}.bam" 2>>"$LOGFILE"
                samtools index "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}.sorted.bam" 2>>"$LOGFILE"
                rm -f "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_R1.sai" "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_R2.sai" "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}.sam" "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}.bam" 2>>"$LOGFILE"

                damageprofiler -i "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}.sorted.bam" -r "$REF_FASTA" -o "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_damageprofiler_unmerged" 2>>"$LOGFILE"
                mapDamage -i "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}.sorted.bam" -r "$REF_FASTA" --folder "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_mapDamage_unmerged" 2>>"$LOGFILE"
            else
                echo "WARNING: Missing output FASTQ files for $GROUP unmerged pairs" | tee -a "$LOGFILE"
            fi
        fi

        # Extraction pour reads merged
        if [[ -f "$MERGED_FILE" ]]; then
            echo "Extracting reads for $GROUP (merged)..." | tee -a "$LOGFILE"
            python3 /home/plstenge/KrakenTools/extract_kraken_reads.py \
                -k "$KRAKEN_FILE" -s "$MERGED_FILE" -t "$TAX_ID" -o "$OUT_MERGED" --fastq-output 2>>"$LOGFILE"

            if [[ -f "$OUT_MERGED" ]]; then
                echo "Aligning with bwa aln/samse (merged)..." | tee -a "$LOGFILE"
                bwa aln -t 4 "$REF_FASTA" "$OUT_MERGED" > "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.sai" 2>>"$LOGFILE"
                bwa samse "$REF_FASTA" "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.sai" "$OUT_MERGED" > "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.sam" 2>>"$LOGFILE"
                samtools view -bS "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.sam" > "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.bam" 2>>"$LOGFILE"
                samtools sort -o "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.sorted.bam" "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.bam" 2>>"$LOGFILE"
                samtools index "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.sorted.bam" 2>>"$LOGFILE"
                rm -f "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.sai" "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.sam" "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.bam" 2>>"$LOGFILE"

                damageprofiler -i "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.sorted.bam" -r "$REF_FASTA" -o "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_damageprofiler_merged" 2>>"$LOGFILE"
                mapDamage -i "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.sorted.bam" -r "$REF_FASTA" --folder "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_mapDamage_merged" 2>>"$LOGFILE"
            else
                echo "WARNING: Missing output FASTQ file for $GROUP merged" | tee -a "$LOGFILE"
            fi
        fi
    done
done

echo "Finished at $(date)" | tee -a "$LOGFILE"
