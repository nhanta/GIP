#!/bin/bash
ref_genome=$1
# Generate a chromosome renaming file
for CHR in {1..23} X ; do 
    echo ${CHR} chr${CHR}
done >> chr_names.txt

# In order to preserve them throughout the protocol, 
#set ID field with unique IDs e.g. in format CHR_POS_REF_ALT
for CHR in {1..22} X; do
    bcftools annotate \
    --set-id '%CHROM\_%POS\_%REF\_%ALT' \
    ../ref_panel/ALL.chr${CHR}_GRCh38.genotypes.20170504.vcf.gz  \
    -Oz -o panel_SNPID_chr${CHR}.vcf.gz
done

# Multiple processing commands piped together
for CHR in {1..22} X; do
    bcftools annotate --rename-chrs chr_names.txt \
        panel_SNPID_chr${CHR}.vcf.gz -Ou | \
    bcftools view -i 'INFO/VT="SNP" | INFO/VT="INDEL"' -Ou | \
    bcftools norm -f ${ref_genome} -d none -Ou | \
    bcftools view -g ^miss -Oz -o 1000GP_chr${CHR}.vcf.gz
done

# Fix the chromosome X ploidy to phased diploid
# Requires a ploidy.txt file containing 
# space-separated CHROM,FROM,TO,SEX,PLOIDY 
echo "chrX 1 156040895 M 2" > ploidy.txt
bcftools +fixploidy \
    1000GP_chrX.vcf.gz -Ov -- -p ploidy.txt | \
    sed 's#0/0#0\|0#g;s#1/1#1\|1#g' | \
bcftools view -Oz -o 1000GP_chr23.vcf.gz

# Duplicate ID removal Remove duplicate IDs
for CHR in {1..23}; do
    bcftools query -f '%ID\n' 1000GP_chr${CHR}.vcf.gz | \
    sort | uniq -d > 1000GP_chr${CHR}.dup_id

    if [[ -s 1000GP_chr${CHR}.dup_id ]]; then
    	bcftools view -e ID=@1000GP_chr${CHR}.dup_id \
    	1000GP_chr${CHR}.vcf.gz \
        -Oz -o 1000GP_filtered_chr${CHR}.vcf.gz
    else 
    	mv 1000GP_chr${CHR}.vcf.gz \
        1000GP_filtered_chr${CHR}.vcf.gz
    fi
done

# Check if the VCF does NOT contain AF in the INFO field,
# and calculate it with bcftools +fill-tags plugin
for CHR in {1..23}; do
    bcftools +fill-tags \
        1000GP_filtered_chr${CHR}.vcf.gz \
        -Oz -o 1000GP_AF_chr${CHR}.vcf.gz -- -t AF
done

# Generate a tab-delimited header
echo -e 'CHR\tSNP\tREF\tALT\tAF' \
    > 1000GP_imputation_all.frq

# Query the required fields from the VCF file 
# and append to the allele frequency file 
for CHR in {1..23}; do
    bcftools query \
    -f '%CHROM\t%CHROM\_%POS\_%REF\_%ALT\t%REF\t%ALT\t%INFO/AF\n' \
    1000GP_AF_chr${CHR}.vcf.gz \
    >> 1000GP_imputation_all.frq
done

bcftools query -l 1000GP_AF_chr22.vcf.gz \
    > 1000GP_sample_IDs.txt