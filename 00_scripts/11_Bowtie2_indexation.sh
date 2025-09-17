#!/bin/bash
#SBATCH --job-name=11_Bowtie2_indexation
#SBATCH --ntasks=1
#SBATCH -p smp
#SBATCH --mem=1000G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL 
#SBATCH --error="/home/plstenge/seda_DNA_Corsican_wreck_HOPS/00_scripts/11_Bowtie2_indexation.err"
#SBATCH --output="/home/plstenge/seda_DNA_Corsican_wreck_HOPS/00_scripts/11_Bowtie2_indexation.out"

# Activation environnement
module load conda/4.12.0
source ~/.bashrc
conda activate metaDMG

#bowtie2-build /storage/biodatabanks/ncbi/NT/current/fasta/All/all.fasta /home/plstenge/refmetagenome

input_fasta="/storage/biodatabanks/ncbi/NT/current/fasta/All/all.fasta"
output_fasta="/home/plstenge/all_NT_hpc2_filtered_for_bowtie2.fasta"

awk '
BEGIN {RS=">"; ORS=""} 
NR>1 {
  header = substr($1, 1, index($1, "\n")-1)
  seq = substr($1, index($1, "\n")+1)
  gsub("\n", "", seq)
  # Remove gaps and Ns
  filtered = seq
  gsub("[-Nn]", "", filtered)
  # If after removing gaps and Ns there is at least one nucleotide left, print sequence
  if(length(filtered) > 0) {
    print ">" header "\n"
    # Print sequence in lines of 60 characters maximum
    for(i=1; i<=length(seq); i+=60) {
      print substr(seq, i, 60) "\n"
    }
  }
}
' $input_fasta > $output_fasta

echo "Filtered FASTA saved to $output_fasta"


# Adapter le chemin selon la disponibilité du scratch
SCRATCHDIR=/home/plstenge/bowtie2_nt_index
mkdir -p $SCRATCHDIR

bowtie2-build --large-index $output_fasta $SCRATCHDIR/refmetagenome
# => Les fichiers $SCRATCHDIR/refmetagenome.*.bt2*
