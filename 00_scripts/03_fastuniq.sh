#!/bin/bash
#SBATCH --job-name=02_adapterremoval
#SBATCH --ntasks=1
#SBATCH -p smp
#SBATCH --mem=250G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL 
#SBATCH --error="/home/plstenge/seda_DNA_Corsican_wreck_HOPS/00_scripts/02_adapterremoval.err"
#SBATCH --output="/home/plstenge/seda_DNA_Corsican_wreck_HOPS/00_scripts/02_adapterremoval.out"

INPUT="/home/plstenge/seda_DNA_Corsican_wreck_HOPS/03_bbduk"
OUTPUT="/home/plstenge/seda_DNA_Corsican_wreck_HOPS/04_fastuniq"

mkdir -p "$OUTPUT"

module load conda/4.12.0
source ~/.bashrc
conda activate fastuniq

cd "$INPUT" || exit 1

# Boucle sur les fichiers R1
for R1 in clean_*_R1.fastq.gz; do
    # En déduire le nom du sample sans suffixe
    base=$(echo "$R1" | sed 's/_R1\.fastq\.gz//')
    R2="${base}_R2.fastq.gz"
    
    # Vérifie que le fichier R2 existe
    if [[ -f "$R2" ]]; then
        echo "Traitement de la paire: $base"
        
        # Crée le fichier .list
        listfile="${OUTPUT}/${base}.list"
        echo -e "${INPUT}/${R1}\n${INPUT}/${R2}" > "$listfile"
        
        # Lancer fastuniq
        fastuniq -i "$listfile" -t q \
            -o "${OUTPUT}/${base}_dedup_R1.fastq" \
            -p "${OUTPUT}/${base}_dedup_R2.fastq"
    else
        echo "ATTENTION: fichier R2 manquant pour $base"
    fi
done

echo "Terminé."
