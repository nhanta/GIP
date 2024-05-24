#!/bin/bash
input_path=$1 
data_name=$2 # Prefix name of data

# Unzip imputed haplotype data
gunzip ${input_path}/output/${data_name}.missing.imputed.hap.gz

# Diplotype to haplotype
bcftools convert --haplegendsample ${input_path}/${data_name}.ori  ${input_path}/${data_name}.ori.vcf
gunzip ${input_path}/${data_name}.ori.hap.gz

echo "GIP"
# Evaluate gip
python gip/evaluation.py \
${input_path}/${data_name}.ori.hap \
${input_path}/${data_name}.missing.hap \
${input_path}/output/${data_name}.missing.imputed.hap

# Diplotype to haplotype
bcftools convert --haplegendsample ${input_path}/output/${data_name}.missing.imputed.beagle ${input_path}/output/${data_name}.missing.imputed.beagle.vcf
gunzip ${input_path}/output/${data_name}.missing.imputed.beagle.hap.gz

echo "BEAGLE5.4"
# Evaluate Beagle5.4
python gip/evaluation.py \
${input_path}/${data_name}.ori.hap \
${input_path}/${data_name}.missing.hap \
${input_path}/output/${data_name}.missing.imputed.beagle.hap