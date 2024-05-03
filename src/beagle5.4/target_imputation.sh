#!/bin/bash
ref_genome_dir=$1 # /work/users/minhnth/gatk/hg38.fasta
input_dir=$2 # ../data/obesity
file_name=$3
output_dir=$4 # ../data/obesity/imputation_output/
ref_panel=$5
genetic_map=$6
num_threads=$7

mkdir ${output_dir}

# Convert multi allelic to bi allelic
bcftools norm \
-m-any \
--check-ref \
-w \
-f ${ref_genome_dir} \
${input_dir}/${file_name} \
-o ${output_dir}/biallelic_snps.vcf

# A comma-separated string of chrs to keep
chrs=$(echo chr{1..22} chrX | tr ' ' ',')

# Keep only those wanted chromosomes
bcftools view -t $chrs ${output_dir}/biallelic_snps.vcf \
    -Oz -o ${output_dir}/biallelic_chrfiltered.vcf.gz

# Align the alleles to the reference genome, 
# and keep only biallelic records
bcftools norm -f ${ref_genome_dir} -c ws \
    ${output_dir}/biallelic_chrfiltered.vcf.gz -Ou | \
bcftools view -m 2 -M 2 \
    -Oz -o ${output_dir}/biallelic_refcorrected.vcf.gz

# Replace the ID column with a CHR_POS_REF_ALT
bcftools annotate \
    --set-id '%CHROM\_%POS\_%REF\_%ALT' \
    ${output_dir}/biallelic_refcorrected.vcf.gz \
    -Oz -o ${output_dir}/biallelic_SNPID.vcf.gz

# Copy the panel sample ID file with 
# a different name to the working directory
cp ${ref_panel}/1000GP_sample_IDs.txt ${output_dir}/duplicate_sample_IDs.txt
cp ${ref_panel}/1000GP_imputation_all.frq ${output_dir}

# Generate a list of chip data sample IDs, 
# keep only duplicates and 
# append to the list of panel sample IDs
bcftools query -l ${output_dir}/biallelic_SNPID.vcf.gz | \
    uniq -d >> ${output_dir}/duplicate_sample_IDs.txt

# Exclude duplicate samples and MAC=0 variants, 
# and update AF value after sample removals.
# Finally, exclude duplicate variants
bcftools view -S ^${output_dir}/duplicate_sample_IDs.txt \
   --force-samples ${output_dir}/biallelic_SNPID.vcf.gz -Ou | \
#bcftools +fill-tags -Ou -- -t AC,AN,AF | \
#bcftools norm -d none \
    #-Oz -o  ${output_dir}/biallelic_noduplicate_variants.vcf.gz
#Don't use bcftools norm for ref_pannel input to addvoid removing all genotypes
bcftools +fill-tags -Ou -- -t AC,AN,AF 
cp ${output_dir}/biallelic_SNPID.vcf.gz ${output_dir}/biallelic_noduplicate_variants.vcf.gz
# Remove low allele count variants 
# (redundant if they are already removed)
# we skip this step for considering rare variants
: << 'COMMENT'
bcftools view -e 'INFO/AC<3 | (INFO/AN-INFO/AC)<3' \
    ${output_dir}/biallelic_noduplicate_variants.vcf.gz \
    -Oz -o ${output_dir}/biallelic_AF.vcf.gz
COMMENT

# Generate a tab-delimited header
echo -e 'CHR\tSNP\tREF\tALT\tAF' \
    > ${output_dir}/biallelic_chip.frq

# Query the required fields from the VCF file 
# and append to the allele frequency file
bcftools query \
    -f '%CHROM\t%ID\t%REF\t%ALT\t%INFO/AF\n' \
    ${output_dir}/biallelic_noduplicate_variants.vcf.gz \
    >> ${output_dir}/biallelic_chip.frq

# Copy and save the given 'plot_AF.R' file and run it with:
Rscript --no-save plot_AF.R \
    ${output_dir}/biallelic_chip.frq \
    ${output_dir}/biallelic \
    ${output_dir}/1000GP_imputation_all.frq \
    0.1 \
    5 \
    ${output_dir}


# Sort the exclusion lists
sort -V ${output_dir}/biallelic_exclude.txt \
    > ${output_dir}/biallelic_exclusion1.txt
sort -V ${output_dir}/biallelic_nonpanel_exclude.txt \
    > ${output_dir}/biallelic_exclusion2.txt

# Exclude the variants showing highly discrepant AF
bcftools view -e ID=@${output_dir}/biallelic_exclusion1.txt  \
    ${output_dir}/biallelic_noduplicate_variants.vcf.gz \
    -Oz -o ${output_dir}/biallelic_exclusion1.vcf.gz

# Exclude the variants not present in the reference panel
bcftools view -e ID=@${output_dir}/biallelic_exclusion2.txt \
    ${output_dir}/biallelic_exclusion1.vcf.gz \
    -Oz -o ${output_dir}/biallelic_for_phasing.vcf.gz

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
        gt=${output_dir}/biallelic_for_phasing.vcf.gz \
        ref=${ref_panel}/1000GP_AF_chr${CHR}.vcf.gz \
        map=${genetic_map}/beagle_chr${CHR}_b38.map \
        out=${output_dir}/biallelic_phased_chr${CHR} impute=false gp=true \
        nthreads=${num_threads} window=500.0 chrom=chr${chrname} \
	    iterations=20 ne=20000 seed=-99999 
	
done

# Re-calculate and add INFO field values for each chromosome
for CHR in {1..23}; do

    # Create index file
    tabix -fp vcf ${output_dir}/biallelic_phased_chr${CHR}.vcf.gz
	
    # Re-calculate allele frequency and 
    # compute Impute2-like INFO score
    bcftools +fill-tags ${output_dir}/biallelic_phased_chr${CHR}.vcf.gz \
        -Ou -- -t AF | \
    bcftools +impute-info \
        -Oz -o ${output_dir}/biallelic_imputed_info_chr${CHR}.vcf.gz
done

# Create imputed list for all chromosomes
inputed_files=${output_dir}'/imputed_list.txt'
chr_names=${output_dir}'/chr_names.txt'

for CHR in {1..23}; do
    if [ -f ${output_dir}/biallelic_imputed_info_chr${CHR}.vcf.gz ]; then
        echo ${output_dir}/biallelic_imputed_info_chr${CHR}.vcf.gz >> $inputed_files
        echo ${CHR} >> ${chr_names}
    else
        echo biallelic_imputed_info_chr${CHR}.vcf.gz 'does not exist.'
    fi
done

# Concat imputed files

bcftools concat -f ${inputed_files} -Oz -o ${output_dir}/biallelic_combination.vcf
rm ${inputed_files}

# Sort SNPs
bcftools sort ${output_dir}/biallelic_combination.vcf -Oz -o ${input_dir}/imputed_biallelic_geno.vcf

: << 'COMMENT'
# List existing chromosomes
chr=()
while IFS= read -r line
    do
        chr+=($(($line)))
    done < ${chr_names}
rm ${chr_names}

# Generate an allele frequency file for plotting for each chrs
for CHR in ${chr}; do
    # Generate a header for the output file
    echo -e 'CHR\tSNP\tREF\tALT\tAF\tINFO\tAF_GROUP' \
        > ${output_dir}/biallelic_group_chr${CHR}.txt

    # Query the required fields and 
    # add frequency group (1, 2 or 3) as the last column
    bcftools query -f \
        '%CHROM\t%CHROM\_%POS\_%REF\_%ALT\t%REF\t%ALT\t%INFO/AF\t%INFO/INFO\t-\n' \
        ${output_dir}/biallelic_imputed_info_chr${CHR}.vcf.gz | \
    # Here $5 refers to AF values, $7 refers to AF group
    awk -v OFS="\t" \
        '{if ($5>=0.05 && $5<=0.95) $7=1; \
            else if(($5>=0.005 && $5<0.05) || \
            ($5<=0.995 && $5>0.95)) $7=2; else $7=3} \
            { print $1, $2, $3, $4, $5, $6, $7 }' \
        >> ${output_dir}/biallelic_group_chr${CHR}.txt
done

# Copy and save the given 'plot_INFO_and_AF_for_imputed_chrs.R' file 
# and run it with:
for CHR in ${chr}; do
    Rscript --no-save  \
        plot_INFO_and_AF_for_imputed_chrs.R \
        ${output_dir}/biallelic_group_chr${CHR}.txt \
        ${output_dir}/biallelic_chr${CHR} \
        ${output_dir}/1000GP_imputation_all.frq \
        ${output_dir}
done

# Combine the plots per chromosome into a single pdf file
convert $(ls ${output_dir}/_chr*.png | sort -V) \
     ${output_dir}.pdf

COMMENT