#!/bin/bash
ref_genome_dir=$1 # /work/users/minhnth/gatk/hg38.fasta
input_dir=$2 # ../data/obesity
file_name=$3 # prefix of filename
output_dir=$4 # ../data/obesity/imputation_output/
ref_panel=$5
genetic_map=$6
num_threads=$7

mkdir -p ${output_dir}

echo "Waiting for Beagle5.4"
log=$(readlink -f ${file_name}.log.txt)
SECONDS=0
# Phasing by Beagle
for CHR in {1..23}; do 
    # Chromosome 23 is coded as 'chrX', 
    # whereas autosomal chrs are numbers alone
    if [ "$CHR" == "23" ]; then
          chrname=X
    else
          chrname=$CHR
    fi

    # Run phasing for each chromosome separately
    beagle -Xmx96g \
        gt=${input_dir}/${file_name}.vcf \
        ref=${ref_panel}/1000GP_AF_chr${CHR}.vcf.gz \
        map=${genetic_map}/beagle_chr${CHR}_b38.map \
        out=${output_dir}/${file_name}.phased.chr${CHR} impute=false gp=true \
        nthreads=${num_threads} window=500.0 chrom=chr${chrname} \
	    iterations=20 ne=20000 seed=-99999 
	
done

# Re-calculate and add INFO field values for each chromosome
for CHR in {1..23}; do

    # Create index file
    tabix -fp vcf ${output_dir}/${file_name}.phased.chr${CHR}.vcf.gz
	
    # Re-calculate allele frequency and 
    # compute Impute2-like INFO score
    bcftools +fill-tags ${output_dir}/${file_name}.phased.chr${CHR}.vcf.gz \
        -Ou -- -t AF | \
    bcftools +impute-info \
        -Oz -o ${output_dir}/${file_name}.imputed_info_chr${CHR}.vcf.gz
done

# Create imputed list for all chromosomes
inputed_files=${output_dir}'/imputed_list.txt'
chr_names=${output_dir}'/chr_names.txt'

for CHR in {1..23}; do
    if [ -f ${output_dir}/${file_name}.imputed_info_chr${CHR}.vcf.gz ]; then
        echo ${output_dir}/${file_name}.imputed_info_chr${CHR}.vcf.gz >> $inputed_files
        echo ${CHR} >> ${chr_names}
    else
        echo ${output_dir}/${file_name}.imputed_info_chr${CHR}.vcf.gz 'does not exist.'
    fi
done

# Concat imputed files

bcftools concat -f ${inputed_files} -Oz -o ${output_dir}/${file_name}.combination.vcf
rm ${inputed_files}

# Sort SNPs
bcftools sort ${output_dir}/${file_name}.combination.vcf -Oz -o ${output_dir}/${file_name}.imputed.beagle.vcf

duration=$SECONDS
echo "$((duration / 60)) minutes and $((duration % 60)) seconds elapsed."
echo "Sample: ${file_name}" >> $log
echo "$((duration / 60)) minutes and $((duration % 60)) seconds elapsed." >> $log