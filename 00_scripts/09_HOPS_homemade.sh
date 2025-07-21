#!/bin/bash
#SBATCH --job-name=07_kraken
#SBATCH --ntasks=1
#SBATCH -p smp
#SBATCH --mem=1000G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL 
#SBATCH --error="/home/plstenge/seda_DNA_Corsican_wreck_HOPS/00_scripts/07_kraken.err"
#SBATCH --output="/home/plstenge/seda_DNA_Corsican_wreck_HOPS/00_scripts/07_kraken.out"

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

# Pour chaque fichier Kraken (.report)
for REPORT in $KRAKEN_DIR/*.report; do
    SAMPLE=$(basename $REPORT .report)
    SAMPLE=${SAMPLE%_*_*_*} # Simplification du nom de l'échantillon

    # Pour chaque groupe taxonomique
    for GROUP in "${!TAXONS[@]}"; do
        TAX_ID="${TAXONS[$GROUP]%:*}"           # Récupère le tax_id avant ':'
        REF_FASTA="${TAXONS[$GROUP]#*:}"        # Récupère le chemin du génome après ':'

        echo "Extracting $GROUP ($TAX_ID) reads from $SAMPLE..."
        python3 extract_kraken_reads.py \
            -k $KRAKEN_DIR/${SAMPLE}.kraken \
            -s $FASTQ_DIR/${SAMPLE}_R1.fastq.gz \
            -s2 $FASTQ_DIR/${SAMPLE}_R2.fastq.gz \
            -t $TAX_ID \
            -o $DAMAGE_DIR/${SAMPLE}_${GROUP}_R1.fastq \
            -o2 $DAMAGE_DIR/${SAMPLE}_${GROUP}_R2.fastq \
            --fastq-output

        if [ -e $DAMAGE_DIR/${SAMPLE}_${GROUP}_R1.fastq ]; then
            echo "Aligning $GROUP reads from $SAMPLE..."
            bwa mem -t 4 $REF_FASTA \
                $DAMAGE_DIR/${SAMPLE}_${GROUP}_R1.fastq \
                $DAMAGE_DIR/${SAMPLE}_${GROUP}_R2.fastq > $DAMAGE_DIR/${SAMPLE}_${GROUP}.sam
            samtools view -bS $DAMAGE_DIR/${SAMPLE}_${GROUP}.sam > $DAMAGE_DIR/${SAMPLE}_${GROUP}.bam
        fi

        # 3. Calculer le profil de dommage
        if [ -e $DAMAGE_DIR/${SAMPLE}_${GROUP}.bam ]; then
            echo "Running damage profile for $GROUP on $SAMPLE..."
            damageprofiler -i $DAMAGE_DIR/${SAMPLE}_${GROUP}.bam \
                  -r $REF_FASTA \
                  -o $DAMAGE_DIR/${SAMPLE}_${GROUP}_damage
        fi
    done
done
done
