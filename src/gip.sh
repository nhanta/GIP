#!/bin/bash
data_path=$1 
data_name=$2 # Prefix name of data
method=$3 # GAIN or GIP
num_al=$4
num_threads=$5

# Create output directory
output=${data_path}/output
mkdir ${output}

if [[ ${num_al} == 2 ]]; then
    # Remove existing files
    rm ${output}/${data_name}.imputed.${method}.hap
    rm ${data_path}/${data_name}.hap
    rm ${data_path}/${data_name:0:-7}ori.hap
    # Diplotype to haplotype
    bcftools convert --haplegendsample ${data_path}/$data_name ${data_path}/${data_name}.vcf
    gunzip -f ${data_path}/${data_name}.hap.gz

    log=$(readlink -f ${data_name}.log.txt)
    SECONDS=0
    # Run imputation 
    python3 gip/main.py --missing_data ${data_path}/${data_name}.hap \
    --method ${method} \
    --num_al ${num_al} \
    --output_data ${output}/${data_name}.imputed.${method}.hap \
    --batch_size 10 --hint_rate 0.9 --alpha 100 \
    --iterations 10000 \
    --num_threads ${num_threads}

    duration=$SECONDS
    echo "$((duration / 60)) minutes and $((duration % 60)) seconds elapsed."
    echo "Sample: ${data_name}" >> $log
    echo "$((duration / 60)) minutes and $((duration % 60)) seconds elapsed." >> $log

    # Export imputed data
    bgzip -f ${output}/${data_name}.imputed.${method}.hap
    cp ${data_path}/${data_name}.samples  ${output}/${data_name}.imputed.${method}.samples
    cp ${data_path}/${data_name}.legend.gz  ${output}/${data_name}.imputed.${method}.legend.gz

    bcftools convert --haplegendsample2vcf ${output}/${data_name}.imputed.${method} > ${output}/${data_name}.imputed.${method}.vcf
else
    # Remove existing files
    rm ${output}/${data_name}.imputed.${method}.hap
    rm ${data_path}/${data_name}.hap
    rm ${data_path}/${data_name:0:-7}ori.hap

    # Generate header
    cat ${data_path}/$data_name.vcf | grep '^##' > ${data_path}/$data_name.header.txt

    log=$(readlink -f ${data_name}.log.txt)
    SECONDS=0
    # Run imputation 
    python3 gip/main.py --missing_data ${data_path}/${data_name}.hap \
    --method ${method} \
    --num_al ${num_al} \
    --output_data ${output}/${data_name}.imputed.${method}.hap \
    --batch_size 10 --hint_rate 0.9 --alpha 100 \
    --iterations 10000 \
    --num_threads ${num_threads}
    
    duration=$SECONDS
    echo "$((duration / 60)) minutes and $((duration % 60)) seconds elapsed."
    echo "Sample: ${data_name}" >> $log
    echo "$((duration / 60)) minutes and $((duration % 60)) seconds elapsed." >> $log
fi


