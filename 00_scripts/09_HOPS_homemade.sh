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

# Boucle sur TOUS les fichiers .kraken (pas .report !)
for KRAKEN_FILE in "$KRAKEN_DIR"/*.kraken; do
    KRAKEN_BASE=$(basename "$KRAKEN_FILE" .kraken)
    # On veut extraire le préfixe commun entre .kraken et .fastq.gz
    # Ex: clean_1121_sed8_rep1_dedup_clumpify_unmerged.kraken → clean_1121_sed8_rep1
    # Le préfixe est tout sauf _dedup_clumpify_unmerged ou _dedup_clumpify_merged
    PREFIX=$(echo "$KRAKEN_BASE" | sed -E 's/_dedup_clumpify_(un)?merged$//')
    R1_FILE="${FASTQ_DIR}/${PREFIX}_dedup_clumpify_fastp_R1.fastq.gz"
    R2_FILE="${FASTQ_DIR}/${PREFIX}_dedup_clumpify_fastp_R2.fastq.gz"
    
    # Pour chaque taxon d'intérêt
    for GROUP in "${!TAXONS[@]}"; do
        TAX_ID="${TAXONS[$GROUP]%:*}"
        REF_FASTA="${TAXONS[$GROUP]#*:}"
        
        OUT_R1="${DAMAGE_DIR}/${PREFIX}_${GROUP}_R1.fastq"
        OUT_R2="${DAMAGE_DIR}/${PREFIX}_${GROUP}_R2.fastq"
        
        echo "======= Extracting $GROUP reads from $KRAKEN_BASE ======="
        python3 /home/plstenge/KrakenTools/extract_kraken_reads.py \
            -k "$KRAKEN_FILE" \
            -s "$R1_FILE" \
            -s2 "$R2_FILE" \
            -t "$TAX_ID" \
            -o "$OUT_R1" \
            -o2 "$OUT_R2" \
            --fastq-output
        
        if [[ -e "$OUT_R1" && -e "$OUT_R2" ]]; then
            echo "======= Aligning $GROUP reads from $PREFIX ======="
            bwa mem -t 4 "$REF_FASTA" "$OUT_R1" "$OUT_R2" > "${DAMAGE_DIR}/${PREFIX}_${GROUP}.sam"
            samtools view -bS "${DAMAGE_DIR}/${PREFIX}_${GROUP}.sam" > "${DAMAGE_DIR}/${PREFIX}_${GROUP}.bam"
            
            if [[ -e "${DAMAGE_DIR}/${PREFIX}_${GROUP}.bam" ]]; then
                echo "======= Running damageprofiler for $GROUP from $PREFIX ======="
                damageprofiler -i "${DAMAGE_DIR}/${PREFIX}_${GROUP}.bam" \
                    -r "$REF_FASTA" \
                    -o "${DAMAGE_DIR}/${PREFIX}_${GROUP}_damage"
            fi
        fi
    done
done
