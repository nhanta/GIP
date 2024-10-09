code=('all' 'AFR' 'AMR' 'EAS' 'EUR' 'SAS')
rate=("0.02" "0.05" "0.1" "0.2" "0.25")
# for 2 alleles
python3 data_preparation.py -p ../../data/HLA \
  -s "${code[@]}" \
  -r "${rate[@]}" \
  -m "igsr-1000 genomes on grch38.tsv" -n "HLA1_chr6"
# for 3 alleles
python3 data_preparation.py -p ../../data/HLA_3_alleles \
  -s "${code[@]}" \
  -r "${rate[@]}" \
  -m "igsr-1000 genomes on grch38.tsv" -n "HLA1_chr6"