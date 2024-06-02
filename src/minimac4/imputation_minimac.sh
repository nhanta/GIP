#!/bin/bash
ref_panel=$1 # /work/users/minhnth/obesity/imputation/ref_panel_for_minimac/1000GP_AF_chr6.msav
input_dir=$2 # ../data/HLA_3_alleles/output/HLA1_chr6.SAS.0.25.missing.imputed.minimac.vcf
file_name=$3 # prefix of filename
output_dir=$4 # ../data/obesity/imputation_output/
region=$5 # chr6:292100-351355
genetic_map=$6
num_threads=$7

bgzip -fc ${input_dir}/${file_name}.vcf > ${input_dir}/${file_name}.vcf.gz
tabix -p vcf ${input_dir}/${file_name}.vcf.gz

echo "Waiting for Minimac4"

echo "" > ${file_name}.log.txt
log=$(readlink -f ${file_name}.log.txt)

SECONDS=0
# Run minimac4
minimac4 ${ref_panel} \
${input_dir}/${file_name}.vcf.gz \
-o ${output_dir}/${file_name}.imputed.minimac.vcf \
-r ${region} \
-m ${genetic_map} \
-f GT \
-t 70 
duration=$SECONDS
echo "$((duration / 60)) minutes and $((duration % 60)) seconds elapsed."

echo "Sample: ${file_name}" >> $log
echo "$((duration / 60)) minutes and $((duration % 60)) seconds elapsed." >> $log

# Get haploid data from minimac4's result
python ../gip/vcf2hap.py \
${input_dir}/${file_name:0:-7}ori.vcf \
${input_dir}/${file_name}.vcf \
${output_dir}/${file_name}.imputed.minimac.vcf \
${output_dir}/${file_name}.imputed.minimac.hap \
minimac