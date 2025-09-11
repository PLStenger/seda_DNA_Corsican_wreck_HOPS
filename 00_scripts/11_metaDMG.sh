#!/bin/bash
#SBATCH --job-name=11_metaDMG
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH -p smp
#SBATCH --mem=500G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error="/home/plstenge/seda_DNA_Corsican_wreck_HOPS/00_scripts/11_metaDMG.err"
#SBATCH --output="/home/plstenge/seda_DNA_Corsican_wreck_HOPS/00_scripts/11_metaDMG.out"

# Activation environnement
module load conda/4.12.0
source ~/.bashrc
conda activate metaDMG

WORKDIR=/home/plstenge/seda_DNA_Corsican_wreck_HOPS/06_fastp
REFIDX=/home/plstenge/refmetagenome
TAXDIR=/home/plstenge/ncbi_tax_dmp
OUTDIR=/home/plstenge/seda_DNA_Corsican_wreck_HOPS/11_metaDMG

mkdir -p "$OUTDIR"


for fqmerged in $WORKDIR/*_merged.fastq*; do
  BASENAME=$(basename "$fqmerged" _dedup_clumpify_fastp_merged.fastq)
  # Alignement
  bowtie2 -x $REFIDX -U $fqmerged -S $WORKDIR/${BASENAME}_merged.sam --very-sensitive-local -p 8
  samtools view -bS $WORKDIR/${BASENAME}_merged.sam | samtools sort -o $WORKDIR/${BASENAME}_merged.sorted.bam
  # LCA metaDMG : estimation dommage par taxon
  metaDMG-cpp lca --threads 8 \
    --bam $WORKDIR/${BASENAME}_merged.sorted.bam \
    --names $TAXDIR/names.dmp \
    --nodes $TAXDIR/nodes.dmp \
    --acc2tax $TAXDIR/nucl_gb.accession2taxid \
    --out $OUTDIR/${BASENAME}_merged
  # Statistiques dégâts
  metaDMG-cpp dfit $WORKDIR/${BASENAME}_merged.bdamage.gz \
    --names $TAXDIR/names.dmp \
    --nodes $TAXDIR/nodes.dmp \
    --out $OUTDIR/${BASENAME}_merged.dfit
done


# Répéter bloc pour les non fusionnés si besoin


# Explication détaillées:

#### # 1. Préparation du répertoire taxonomy NCBI (à faire une fois)
#### 
#### mkdir ncbi_tax_dmp
#### cd ncbi_tax_dmp
#### wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip
#### unzip new_taxdump.zip
#### wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
#### gunzip nucl_gb.accession2taxid.gz
#
#### # 2. Alignement des reads sur les génomes cibles
#### # Exemple avec Bowtie2 pour un génome metagenome indexé :
#### # Indexage du génome cible (si pas déjà fait)
#### bowtie2-build /storage/biodatabanks/ncbi/NT/current/fasta/All/all.fasta /home/plstenge/refmetagenome
#
## Pour des reads fusionnés :
#bowtie2 -x refmetagenome -U clean_1120_sed6_rep3_dedup_clumpify_fastp_merged.fastq \
#        -S sample_merged.sam --very-sensitive-local -p 8
#
## Pour des reads non fusionnés (paire) :
#bowtie2 -x refmetagenome -1 clean_1120_sed6_rep3_dedup_clumpify_fastp_R1.fastq \
#       -2 clean_1120_sed6_rep3_dedup_clumpify_fastp_R2.fastq \
#       -S sample_paired.sam --very-sensitive-local -p 8
#
## Convertir ensuite les fichiers SAM en BAM triés :
#
#samtools view -bS sample_merged.sam | samtools sort -o sample_merged.sorted.bam
#samtools view -bS sample_paired.sam | samtools sort -o sample_paired.sorted.bam
#
## 3. metaDMG LCA : assignation taxonomique + matrice de dommages
#
## À adapter pour chaque fichier BAM obtenu :
#
#metaDMG-cpp lca \
#  --threads 8 \
#  --bam sample_merged.sorted.bam \
#  --names ncbi_tax_dmp/names.dmp \
#  --nodes ncbi_tax_dmp/nodes.dmp \
#  --acc2tax ncbi_tax_dmp/nucl_gb.accession2taxid \
#  --sim_score_low 0.95 --sim_score_high 1.0 \
#  --out sample_merged
#
## 4. Estimation du dommage (statistiques de déamination)
#
#metaDMG-cpp dfit sample_merged.bdamage.gz \
#  --names ncbi_tax_dmp/names.dmp \
#  --nodes ncbi_tax_dmp/nodes.dmp \
#  --showfits 0 \
#  --out sample_merged.dfit
#
## Pour des résultats avec intervalles de confiance :
#
#metaDMG-cpp dfit sample_merged.bdamage.gz --nbootstrap 100 --showfits 2 --out sample_merged.dfit
#
## 5. Visualiser ou exploiter les résultats
#
## Pour obtenir la matrice de substitutions ou les détails par taxid/espèce/clade :
#
#metaDMG-cpp print sample_merged.bdamage.gz -names ncbi_tax_dmp/names.dmp
