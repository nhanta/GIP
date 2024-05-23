#!/bin/bash
data_path=$1 
data_name=$2 # Prefix name of data
method=$3 
num_al=$4
num_threads=$5

# Create output directory
output=${data_path}/output
mkdir ${output}

# Diplotype to haplotype
bcftools convert --haplegendsample ${data_path}/$data_name ${data_path}/${data_name}.vcf
gunzip ${data_path}/${data_name}.hap.gz

# Run imputation 
python3 gip/main.py --missing_data ${data_path}/${data_name}.hap \
--method ${method} \
--num_al ${num_al} \
--output_data ${output}/${data_name}.imputed.hap \
--batch_size 10 --hint_rate 0.9 --alpha 100 \
--iterations 10000 \
--num_threads ${num_threads}

# Export imputed data
bgzip ${output}/${data_name}.imputed.hap
cp ${data_path}/${data_name}.samples  ${output}/${data_name}.imputed.samples
cp ${data_path}/${data_name}.legend.gz  ${output}/${data_name}.imputed.legend.gz

bcftools convert --haplegendsample2vcf ${output}/${data_name}.imputed > ${output}/${data_name}.imputed.gip.vcf
gunzip ${output}/${data_name}.imputed.hap.gz