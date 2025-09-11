#!/bin/bash
#SBATCH --job-name=10_table_assignation
#SBATCH --ntasks=1
#SBATCH -p smp
#SBATCH --mem=1000G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL 
#SBATCH --error="/home/plstenge/seda_DNA_Corsican_wreck_HOPS/00_scripts/10_table_assignation.err"
#SBATCH --output="/home/plstenge/seda_DNA_Corsican_wreck_HOPS/00_scripts/10_table_assignation.out"

cd /home/plstenge/seda_DNA_Corsican_wreck_HOPS/07_kraken2/

#git clone https://github.com/jenniferlu717/KrakenTools.git

python3 /home/plstenge/seda_DNA_Corsican_wreck_HOPS/07_kraken2/KrakenTools/kreport2mpa.py -r clean_1121_sed8_rep1_dedup_clumpify_unmerged.report -o clean_1121_sed8_rep1_dedup_clumpify_unmerged.mpa
python3 /home/plstenge/seda_DNA_Corsican_wreck_HOPS/07_kraken2/KrakenTools/kreport2mpa.py -r clean_1122_sed8_rep2_dedup_clumpify_unmerged.report -o clean_1122_sed8_rep2_dedup_clumpify_unmerged.mpa
python3 /home/plstenge/seda_DNA_Corsican_wreck_HOPS/07_kraken2/KrakenTools/kreport2mpa.py -r clean_1129_sed6_rep1_dedup_clumpify_merged.report -o clean_1129_sed6_rep1_dedup_clumpify_merged.mpa
python3 /home/plstenge/seda_DNA_Corsican_wreck_HOPS/07_kraken2/KrakenTools/kreport2mpa.py -r clean_1129_sed6_rep1_dedup_clumpify_unmerged.report -o clean_1129_sed6_rep1_dedup_clumpify_unmerged.mpa
python3 /home/plstenge/seda_DNA_Corsican_wreck_HOPS/07_kraken2/KrakenTools/kreport2mpa.py -r clean_1130_sed6_rep2_dedup_clumpify_merged.report -o clean_1130_sed6_rep2_dedup_clumpify_merged.mpa
python3 /home/plstenge/seda_DNA_Corsican_wreck_HOPS/07_kraken2/KrakenTools/kreport2mpa.py -r clean_1130_sed6_rep2_dedup_clumpify_unmerged.report -o clean_1130_sed6_rep2_dedup_clumpify_unmerged.mpa
python3 /home/plstenge/seda_DNA_Corsican_wreck_HOPS/07_kraken2/KrakenTools/kreport2mpa.py -r clean_1131_sed8_rep3_dedup_clumpify_merged.report -o clean_1131_sed8_rep3_dedup_clumpify_merged.mpa
python3 /home/plstenge/seda_DNA_Corsican_wreck_HOPS/07_kraken2/KrakenTools/kreport2mpa.py -r clean_1131_sed8_rep3_dedup_clumpify_unmerged.report -o clean_1131_sed8_rep3_dedup_clumpify_unmerged.mpa
python3 /home/plstenge/seda_DNA_Corsican_wreck_HOPS/07_kraken2/KrakenTools/kreport2mpa.py -r clean_NTC_sed_dedup_clumpify_merged.report -o clean_NTC_sed_dedup_clumpify_merged.mpa
python3 /home/plstenge/seda_DNA_Corsican_wreck_HOPS/07_kraken2/KrakenTools/kreport2mpa.py -r clean_NTC_sed_dedup_clumpify_unmerged.report -o clean_NTC_sed_dedup_clumpify_unmerged.mpa

python3 /home/plstenge/seda_DNA_Corsican_wreck_HOPS/07_kraken2/KrakenTools/combine_mpa.py -i clean_1121_sed8_rep1_dedup_clumpify_unmerged.mpa clean_1122_sed8_rep2_dedup_clumpify_unmerged.mpa clean_1129_sed6_rep1_dedup_clumpify_merged.mpa clean_1129_sed6_rep1_dedup_clumpify_unmerged.mpa clean_1130_sed6_rep2_dedup_clumpify_merged.mpa clean_1130_sed6_rep2_dedup_clumpify_unmerged.mpa clean_1131_sed8_rep3_dedup_clumpify_merged.mpa clean_1131_sed8_rep3_dedup_clumpify_unmerged.mpa clean_NTC_sed_dedup_clumpify_merged.mpa clean_NTC_sed_dedup_clumpify_unmerged.mpa -o combined_mpa.tsv



