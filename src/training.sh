#!/bin/bash
ref_genome=$1
ref_panel=$2 # /work/users/minhnth/obesity/imputation/ref_panel_for_minimac/1000GP_AF_chr6.msav
input_dir=$3 # ../data/HLA_3_alleles/output/HLA1_chr6.SAS.0.25.missing.imputed.minimac.vcf
prefix_name=$4
region=$5 # chr6:292100-351355
genetic_map=$6
num_threads=$7
num_als=$8

for rate in 0.02 0.05 0.1 0.2 0.25; do
    for p in 'all' 'AFR' 'AMR' 'EAS' 'EUR' 'SAS'; do
        file_name=${prefix_name}.${p}.${rate}.missing
: << 'COMMENT'
        echo "Run GIP"
        bash gip.sh \
        ${input_dir} \
        ${file_name} \
        gip \
        ${num_als} \
        ${num_threads}
COMMENT
#: << 'COMMENT'
        echo "Run GAIN"
        bash gip.sh \
        ${input_dir} \
        ${file_name} \
        gain \
        ${num_als} \
        ${num_threads}
: << 'COMMENT'
        echo "Run Beagle5.4"
        bash beagle5.4/target_imputation.sh \
        ${ref_genome} \
        ${input_dir} \
        ${file_name} \
        ${input_dir}/output/ \
        ${ref_panel} \
        ${genetic_map} \
        ${num_threads}

        echo "Run Minimac4"
        bash minimac4/imputation_minimac.sh \
        ${ref_panel:0:-6}msav \
        ${input_dir} \
        ${file_name} \
        ${input_dir}/output/ \
        ${region} \
        ${genetic_map} \
        ${num_threads}

        echo "Run Impute5"

        bash impute5/impute.sh \
        ${ref_panel} \
        ${input_dir} \
        ${file_name} \
        ${input_dir}/output/ \
        ${region} \
        ${genetic_map:0:-19}chr6.b38.gmap.gz \
        ${num_threads} \
        impute5_v1.2.0/impute5_v1.2.0_static 
COMMENT
    done
done
