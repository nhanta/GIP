# INPUT:    Post-QC chip genotype data in VCFv4.2 format
#           chrX genotypes as diploid genotypes
# example input: DATASET.vcf.gz

# OUTPUT:   DATASET_imputed_info_chr#.vcf.gz (where # is chromosome number)
#           DATASET_postimputation_summary_plots.pdf

# 1.3.2 Minimum quality control
echo "STEP 1: Minimum quality control" 
# Generate a chromosome renaming file
for CHR in {1..23} X ; do 
    echo ${CHR} chr${CHR}
done >> chr_names.txt
# Multiple processing commands piped together
for CHR in {1..22} X; 
do
    # Rename chromosome 
    bcftools annotate --rename-chrs chr_names.txt \
        ALL.chr${CHR}_GRCh38.genotypes.20170504.vcf.gz -Ou | \
    # Remove the rare variants (remove by AC threshold)
    bcftools view -e 'INFO/AC<3 | INFO/AN-INFO/AC<3' -Ou | \
    # Split multiallelic sites to biallelic
    bcftools norm -m -any -Ou | \
    # Keep only SNPs and INDELs
    bcftools view -i 'INFO/VT="SNP" | INFO/VT="INDEL"' -Ou | \
    # Align the variants to reference genome
    bcftools norm -f Homo_sapiens_assembly38.fasta -d none -Ou | \
    # Remove multiallelic records
    bcftools view -m 2 -M 2 -Ou | \
    # Remove missing data
    bcftools view -g ^miss -Oz -o "1000GP_chr${CHR}.vcf.gz"
done
# If have multiallelic. Set ID field with unique IDs e.g. in format CHR_POS_REF_ALT
for CHR in {1..23}; do
    bcftools annotate \
    --set-id '%CHROM\_%POS\_%REF\_%ALT' \
    panel_chr${CHR}.vcf.gz  \
    -Oz -o panel_SNPID_chr${CHR}.vcf.gz
done
# -------------------------------------------------------------------
# 1.3.3 Convert haploid genotypes to homozygous diploids
echo "STEP 2: Convert haploid genotypes to homozygous diploids" 
# Fix the chromosome X ploidy to phased diploid
# Requires a ploidy.txt file containing 
# space-separated CHROM,FROM,TO,SEX,PLOIDY 
echo "chrX 1 156040895 M 2" > ploidy.txt
bcftools +fixploidy \
    1000GP_chrX.vcf.gz -Ov -- -p ploidy.txt | \
    sed 's#0/0#0\|0#g;s#1/1#1\|1#g' | \
bcftools view -Oz -o "1000GP_chr23.vcf.gz"
# ------------------------------------------------------------------
# 1.3.4 Duplicate ID removal
echo "STEP 3: Remove duplicate ID" 
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
# -------------------------------------------------------------------
# 1.3.5 Reference panel allele frequencies
echo "STEP 4: Reference panel allele frequencies"
# Check if the VCF does NOT contain AF in the INFO field,
# and calculate it with bcftools +fill-tags plugin
for CHR in {1..23}; do
    bcftools +fill-tags \
        1000GP_filtered_chr${CHR}.vcf.gz \
        -Oz -o 1000GP_AF_chr${CHR}.vcf.gz -- -t AF
done

# # Generate a tab-delimited header
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

# 1.3.6 Create binary reference panel files\
echo "STEP 5: Create binary reference panel files"
# wget https://faculty.washington.edu/browning/beagle/bref.27Jan18.7e1.jar
# Convert each file to bref format
for CHR in {1..23}; do
    java -jar /path/to/bref.27Jan18.7e1.jar \
        1000GP_AF_chr${CHR}.vcf.gz
done

# 1.3.7 Generate a list of the reference panel sample IDs
echo "STEP 6: Generate a list of the reference panel sample IDs"
bcftools query -l 1000GP_AF_chr22.vcf.gz \
    > "1000GP_sample_IDs.txt"

echo "Done"