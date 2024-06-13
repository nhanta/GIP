#!/bin/bash
input_dir=$1 
prefix_name=$2
num_als=$3

for rate in 0.02 0.05 0.1 0.2 0.25; do
    for p in 'all' 'AFR' 'AMR' 'EAS' 'EUR' 'SAS'; do
        file_name=${prefix_name}.${p}.${rate}
        # Evaluation for 2 alleles
        bash evaluation.sh \
        ${input_dir} \
        ${file_name} \
        ${num_als}
    done
done
