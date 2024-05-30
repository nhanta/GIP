#!/bin/bash
input_path=$1 
data_name=$2 # Prefix name of data
num_al=$3

echo "GIP"
if [[ ${num_al} == 2 ]]; then
    # Unzip imputed haplotype data
    gunzip ${input_path}/output/${data_name}.missing.imputed.gip.hap.gz

    # Diplotype to haplotype
    bcftools convert --haplegendsample ${input_path}/${data_name}.ori  ${input_path}/${data_name}.ori.vcf
    gunzip ${input_path}/${data_name}.ori.hap.gz
fi

# Evaluate gip
python gip/evaluation.py \
${input_path}/${data_name}.ori.hap \
${input_path}/${data_name}.missing.hap \
${input_path}/output/${data_name}.missing.imputed.gip.hap \
${num_al}

echo "BEAGLE5.4"
if [[ ${num_al} == 2 ]]; then
    # Diplotype to haplotype
    bcftools convert --haplegendsample ${input_path}/output/${data_name}.missing.imputed.beagle ${input_path}/output/${data_name}.missing.imputed.beagle.vcf
    gunzip ${input_path}/output/${data_name}.missing.imputed.beagle.hap.gz
fi
# Evaluate Beagle5.4
python gip/evaluation.py \
${input_path}/${data_name}.ori.hap \
${input_path}/${data_name}.missing.hap \
${input_path}/output/${data_name}.missing.imputed.beagle.hap \
${num_al}

echo "MINIMAC4"
python gip/evaluation.py \
${input_path}/${data_name}.ori.hap \
${input_path}/${data_name}.missing.hap \
${input_path}/output/${data_name}.missing.imputed.minimac.hap \
${num_al}