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


function log {
    echo "$(date) -- $1" | tee -a "$LOGFILE"
}

for KRAKEN_FILE in "$KRAKEN_DIR"/*.kraken; do
    KRAKEN_BASE=$(basename "$KRAKEN_FILE" .kraken)
    log "==== Processing file: $KRAKEN_FILE ===="

    # On garde tout le nom de base pour éviter tout problème
    log "Base name for FASTQ search: $KRAKEN_BASE"
    R1_FILE="$FASTQ_DIR/${KRAKEN_BASE}_R1.fastq.gz"
    R2_FILE="$FASTQ_DIR/${KRAKEN_BASE}_R2.fastq.gz"
    MERGED_FILE="$FASTQ_DIR/${KRAKEN_BASE}_merged.fastq.gz"

    log "Looking for R1, R2, merged FASTQ:"
    log "  R1: $R1_FILE"
    log "  R2: $R2_FILE"
    log "  merged: $MERGED_FILE"

    # Pour chaque taxon d'intérêt
    for GROUP in "${!TAXONS[@]}"; do
        TAX_ID="${TAXONS[$GROUP]%:*}"
        REF_FASTA="${TAXONS[$GROUP]#*:}"
        log "--- Processing taxon: $GROUP (taxid: $TAX_ID) ---"

        OUT_R1="$DAMAGE_DIR/${KRAKEN_BASE}_${GROUP}_R1.fastq"
        OUT_R2="$DAMAGE_DIR/${KRAKEN_BASE}_${GROUP}_R2.fastq"
        OUT_MERGED="$DAMAGE_DIR/${KRAKEN_BASE}_${GROUP}_merged.fastq"

        # Extraction pour les paires (unmerged)
        if [[ -f "$R1_FILE" && -f "$R2_FILE" ]]; then
            log "Extracting reads for $GROUP from $KRAKEN_BASE (unmerged pairs)..."
            python3 /home/plstenge/KrakenTools/extract_kraken_reads.py \
                -k "$KRAKEN_FILE" \
                -s "$R1_FILE" \
                -s2 "$R2_FILE" \
                -t "$TAX_ID" \
                -o "$OUT_R1" \
                -o2 "$OUT_R2" \
                --fastq-output

            if [[ -f "$OUT_R1" && -f "$OUT_R2" ]]; then
                log "Read extraction OK. Files present: $OUT_R1, $OUT_R2"
                log "Aligning for $GROUP (unmerged pairs) using bwa aln/sampe..."
                bwa aln -t 4 "$REF_FASTA" "$OUT_R1" > "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_R1.sai" 2>>"$LOGFILE"
                bwa aln -t 4 "$REF_FASTA" "$OUT_R2" > "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_R2.sai" 2>>"$LOGFILE"
                bwa sampe "$REF_FASTA" \
                    "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_R1.sai" \
                    "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_R2.sai" \
                    "$OUT_R1" \
                    "$OUT_R2" > "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}.sam" 2>>"$LOGFILE"
                samtools view -bS "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}.sam" > "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}.bam" 2>>"$LOGFILE"
                samtools sort -o "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}.sorted.bam" "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}.bam" 2>>"$LOGFILE"
                samtools index "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}.sorted.bam" 2>>"$LOGFILE"
                rm "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_R1.sai" "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_R2.sai" "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}.sam" "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}.bam" 2>>"$LOGFILE"
                
                log "Alignment OK. BAM: ${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}.sorted.bam"

                if [[ -f "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}.sorted.bam" ]]; then
                    log "Running damageprofiler for $GROUP (unmerged pairs)..."
                    damageprofiler -i "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}.sorted.bam" \
                        -r "$REF_FASTA" \
                        -o "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_damageprofiler_unmerged" 2>>"$LOGFILE"
                    log "Running mapDamage for $GROUP (unmerged pairs)..."
                    mapDamage -i "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}.sorted.bam" \
                        -r "$REF_FASTA" \
                        --folder "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_mapDamage_unmerged" 2>>"$LOGFILE"
                else
                    log "FATAL: BAM file not generated."
                fi
            else
                log "WARNING: Extraction failed for ${KRAKEN_BASE}_${GROUP} (unmerged pairs). Files missing?"
            fi
        fi

        # Extraction pour les merged reads (single-end)
        if [[ -f "$MERGED_FILE" ]]; then
            log "Extracting reads for $GROUP from $KRAKEN_BASE (merged)..."
            python3 /home/plstenge/KrakenTools/extract_kraken_reads.py \
                -k "$KRAKEN_FILE" \
                -s "$MERGED_FILE" \
                -t "$TAX_ID" \
                -o "$OUT_MERGED" \
                --fastq-output

            if [[ -f "$OUT_MERGED" ]]; then
                log "Read extraction OK (merged). File present: $OUT_MERGED"
                log "Aligning for $GROUP (merged) using bwa aln/samse..."
                bwa aln -t 4 "$REF_FASTA" "$OUT_MERGED" > "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.sai" 2>>"$LOGFILE"
                bwa samse "$REF_FASTA" \
                    "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.sai" \
                    "$OUT_MERGED" > "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.sam" 2>>"$LOGFILE"
                samtools view -bS "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.sam" > "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.bam" 2>>"$LOGFILE"
                samtools sort -o "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.sorted.bam" "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.bam" 2>>"$LOGFILE"
                samtools index "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.sorted.bam" 2>>"$LOGFILE"
                rm "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.sai" "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.sam" "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.bam" 2>>"$LOGFILE"
                
                log "Alignment OK (merged). BAM: ${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.sorted.bam"

                if [[ -f "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.sorted.bam" ]]; then
                    log "Running damageprofiler for $GROUP (merged)..."
                    damageprofiler -i "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.sorted.bam" \
                        -r "$REF_FASTA" \
                        -o "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_damageprofiler_merged" 2>>"$LOGFILE"
                    log "Running mapDamage for $GROUP (merged)..."
                    mapDamage -i "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.sorted.bam" \
                        -r "$REF_FASTA" \
                        --folder "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_mapDamage_merged" 2>>"$LOGFILE"
                else
                    log "FATAL: BAM file not generated (merged)."
                fi
            else
                log "WARNING: Extraction failed for ${KRAKEN_BASE}_${GROUP} (merged). File missing?"
            fi
        fi
    done
done

log "Script finished at $(date)"
