#!/bin/bash
data_path=$1 
data_name=$2 # Prefix name of data
method=$3 
num_al=$4
num_threads=$5

# Create output directory
output=${data_path}/output
mkdir ${output}

# Generate header
cat ${data_path}/$data_name.vcf | grep '^##' > ${data_path}/$data_name.header.txt

# Run imputation 
python3 gip/main.py --missing_data ${data_path}/${data_name}.hap \
--method ${method} \
--num_al ${num_al} \
--output_data ${output}/${data_name}.imputed.hap \
--batch_size 10 --hint_rate 0.9 --alpha 100 \
--iterations 10000 \
--num_threads ${num_threads}
