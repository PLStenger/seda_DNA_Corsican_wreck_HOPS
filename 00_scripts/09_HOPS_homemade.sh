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

for KRAKEN_FILE in "$KRAKEN_DIR"/*.kraken; do
    KRAKEN_BASE=$(basename "$KRAKEN_FILE" .kraken)
    # On NE SIMPLIFIE PAS le nom, car on traite merged et unmerged séparément

    # On cherche les FASTQ correspondants au même nom (attention: ta nomenclature, il faut adapter!)
    # Ex: clean_1121_sed8_rep1_dedup_clumpify_unmerged.kraken → utilise le même nom pour les FASTQ
    R1_FILE="${FASTQ_DIR}/${KRAKEN_BASE}_R1.fastq.gz"
    R2_FILE="${FASTQ_DIR}/${KRAKEN_BASE}_R2.fastq.gz"
    MERGED_FILE="${FASTQ_DIR}/${KRAKEN_BASE}_merged.fastq.gz"  # Si tu as un merged FASTQ

    for GROUP in "${!TAXONS[@]}"; do
        TAX_ID="${TAXONS[$GROUP]%:*}"
        REF_FASTA="${TAXONS[$GROUP]#*:}"

        OUT_R1="${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_R1.fastq"
        OUT_R2="${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_R2.fastq"
        OUT_MERGED="${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.fastq"

        echo "======= Extracting $GROUP reads from $KRAKEN_BASE ======="
        # Extraction pour les paires
        if [[ -f "$R1_FILE" && -f "$R2_FILE" ]]; then
            python3 /home/plstenge/KrakenTools/extract_kraken_reads.py \
                -k "$KRAKEN_FILE" \
                -s "$R1_FILE" \
                -s2 "$R2_FILE" \
                -t "$TAX_ID" \
                -o "$OUT_R1" \
                -o2 "$OUT_R2" \
                --fastq-output

            # ALIGNEMENT BWA ALN (ADN ancien)
            if [[ -f "$OUT_R1" && -f "$OUT_R2" ]]; then
                echo "======= Aligning $GROUP reads (unmerged) from $KRAKEN_BASE (bwa aln/sampe) ======="
                bwa aln -t 4 "$REF_FASTA" "$OUT_R1" > "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_R1.sai"
                bwa aln -t 4 "$REF_FASTA" "$OUT_R2" > "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_R2.sai"
                bwa sampe "$REF_FASTA" \
                    "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_R1.sai" \
                    "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_R2.sai" \
                    "$OUT_R1" \
                    "$OUT_R2" > "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}.sam"
                samtools view -bS "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}.sam" > "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}.bam"
                samtools sort -o "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}.sorted.bam" "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}.bam"
                samtools index "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}.sorted.bam"
                rm "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_R1.sai" "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_R2.sai"

                # DAMAGE PROFILING
                if [[ -f "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}.sorted.bam" ]]; then
                    echo "======= Running damageprofiler for $GROUP (unmerged) ======="
                    damageprofiler -i "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}.sorted.bam" \
                        -r "$REF_FASTA" \
                        -o "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_damageprofiler_unmerged"

                    echo "======= Running mapDamage for $GROUP (unmerged) ======="
                    mapDamage -i "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}.sorted.bam" \
                        -r "$REF_FASTA" \
                        --folder "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_mapDamage_unmerged"
                fi
            fi
        fi

        # Extraction pour les merged
        if [[ -f "$MERGED_FILE" ]]; then
            python3 /home/plstenge/KrakenTools/extract_kraken_reads.py \
                -k "$KRAKEN_FILE" \
                -s "$MERGED_FILE" \
                -t "$TAX_ID" \
                -o "$OUT_MERGED" \
                --fastq-output

            # ALIGNEMENT BWA ALN sur merged (single-end)
            if [[ -f "$OUT_MERGED" ]]; then
                echo "======= Aligning $GROUP reads (merged) from $KRAKEN_BASE (bwa aln/samse) ======="
                bwa aln -t 4 "$REF_FASTA" "$OUT_MERGED" > "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.sai"
                bwa samse "$REF_FASTA" \
                    "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.sai" \
                    "$OUT_MERGED" > "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.sam"
                samtools view -bS "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.sam" > "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.bam"
                samtools sort -o "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.sorted.bam" "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.bam"
                samtools index "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.sorted.bam"
                rm "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.sai"

                # DAMAGE PROFILING sur merged
                if [[ -f "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.sorted.bam" ]]; then
                    echo "======= Running damageprofiler for $GROUP (merged) ======="
                    damageprofiler -i "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.sorted.bam" \
                        -r "$REF_FASTA" \
                        -o "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_damageprofiler_merged"

                    echo "======= Running mapDamage for $GROUP (merged) ======="
                    mapDamage -i "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.sorted.bam" \
                        -r "$REF_FASTA" \
                        --folder "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_mapDamage_merged"
                fi
            fi
        fi
    done
done
