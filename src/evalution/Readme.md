# Evaluation Script Documentation

The script `evaluation.py` evaluates the accuracy of imputed data against original data. Below is a description of each argument:

- `-o`, `--ori_data`: Specifies the path to the original data file (hap format). This argument is required.
    - **Type**: `str`
    - **Example**: `-o /path/to/original_data.hap`

- `-m`, `--missing_data`: Specifies the path to the missing data file (hap format). This argument is required.
    - **Type**: `str`
    - **Example**: `-m /path/to/missing_data.hap`

- `-i`, `--imputed_data`: Specifies the path to the imputed data file (hap format). This argument is required.
    - **Type**: `str`
    - **Example**: `-i /path/to/imputed_data.hap`

- `-p`, `--path`: Specifies the path to the result file where the evaluation results will be appended. This argument is required.
    - **Type**: `int`
    - **Example**: `-p /path/to/result.txt`

- `-m`, `--method`: Specifies the method used for imputation. This argument is required.
    - **Type**: `str`
    - **Example**: `-m Beagle`

You can convert VCF file to HAP file by using `vcf2hap.sh` script in this directory.