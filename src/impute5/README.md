# GIP
GIP: An Advanced Genotype Imputation Method Based on GANs

## Tool Installations
Require conda

```bash
conda env create -f enviroment.yaml
```
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