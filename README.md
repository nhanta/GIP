# GIP
GIP: An Advanced Genotype Imputation Method Based on GANs

## Tool Installations
### Installing GATK 4.2.0.0
```
wget https://github.com/broadinstitute/gatk/releases/download/4.2.0.0/gatk-4.2.0.0.zip
unzip gatk-4.2.0.0.zip -d /path to your directory/
export PATH=$PATH:/path to your directory/
```
### Installing Samtools
```
wget https://sourceforge.net/projects/samtools/files/samtools/
tar xvjf samtools-1.1.tar.bz2
cd samtools-1.1
make
export PATH=$PATH:/path to your directory/
```
### Installing Other Tools
- `conda install -c bioconda bcftools`
- `conda install -c bioconda plink`
- `conda install -c bioconda sra-tools`
- `conda install -c bioconda fastqc`
- `conda install -c bioconda trimmomatic`
- `conda install -c bioconda bwa`
- `conda install -c bioconda tabix`
- `conda install vcftools`
- `pip install bed-reader`
- `pip install captum`
## Data Preparation 
Single-end read data includes [139 obesity samples](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP139885&o=acc_s%3Aa) from Ion Torrent. We run [bash script files](https://github.com/nhanta/Advanced_Methods_for_Disease_Risk_Prediction/tree/main/gatk) to download data, and extract the data to fastq.
- Download pileup data by sratools: `prefetch --option-file SRR_Acc_List.txt`
- Decompress sra to fastq: `sh sra_to_fastq.sh`
- Fix fastq files: `sh fix_fastq.sh`
## Implementing GATK 
In the study, this part belongs to the **initial data preprocessing** and **variant calling**:

`sh gatk_obesity.sh [gatk directory] [input directory] [output directory] [number of threads]`

Folders can be arranged as follows:
```
GATK
|-- Work            
|   |-- GATK directory                   Include GRCh38 and related data.
|   `-- GATK 4.2.0.0                     GATK version 4.2.0.0.                 
|   `-- snpEff/Funcotator                To annotate variants.
|-- Input directory                      Include 139 samples, type of fastq.
|-- Output directory                     To generate result files.
```
## Subsequent data preprocessing
Data after implementing GATK is [transformed to PLINK](https://github.com/nhanta/Advanced_Methods_for_Disease_Risk_Prediction/blob/main/obs_preprocessing.ipynb). Then, we use PLINK to preprocess data following the PRS criteria:
```
plink \
    --bfile obs_ngs \
    --maf 0.05 \
    --hwe 1e-6 \
    --geno 0.1 \
    --write-snplist \
    --make-bed \
    --out obs_ngs.QC
```
Next, we convert PLINK file to one vcf file for training model:

`plink --bfile [filename prefix] --recode vcf --out [VCF prefix]`

## Prediction
We use classical machine learning models, such as [Logistic Regression, Decision Tree, and SVM](https://github.com/nhanta/Advanced_Methods_for_Disease_Risk_Prediction/blob/main/lr_dt_svm.py) to select importance features. Besides, we use a [Neural-Network Feature Selection](https://github.com/nhanta/Advanced_Methods_for_Disease_Risk_Prediction/blob/main/lr_dt_svm.py) for that. 


# Generate missing data
```
python data_preparation.py \
HLA1_chr6 
```

# Imputation by beagle
```
bash target_imputation.sh \
/work/users/minhnth/gatk/hg38.fasta \
../../data/HLA/ \
HLA1_chr6.SAS.0.25.missing \
../../data/HLA/output/ \
/work/users/minhnth/obesity/imputation/ref_panel \
/work/users/minhnth/obesity/imputation/genetic_map_for_imputation \
70
```
# Imputation by GIP
```
bash gip.sh \
../data/HLA \
HLA1_chr6.SAS.0.25.missing \
gip \
2 \
70 
```
# Evaluation for GIP and Beagle 5.4
```
bash evaluation.sh \
../data/HLA/ \
HLA1_chr6.SAS.0.25 
```