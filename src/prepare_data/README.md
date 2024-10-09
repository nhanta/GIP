# Data Preparation
The script `data_preparation.py` will create random data from vcf files. Below is a description of each argument:

- `-p`, `--data_path`: Specifies the path to the directory containing the VCF files and metadata. This argument is required.
    - **Type**: `str`
    - **Example**: `-p /path/to/data`

- `-s`, `--superpopulation_code`: Specifies the super population codes to be selected. This argument is optional and defaults to `['all']`.
    - **Type**: `str`
    - **Default**: `['all']`
    - **Example**: `-s AFR AMR`

- `-r`, `--rate`: Specifies the rate of missing data. This argument is optional and defaults to `['0.25']`.
    - **Type**: `str`
    - **Default**: `['0.25']`
    - **Example**: `-r 0.02 0.05`

- `-m`, `--meta`: Specifies the metadata file. This argument is required.
    - **Type**: `str`
    - **Example**: `-m metadata.tsv`

- `-n`, `--data_name`: Specifies the prefix of the data name. This argument is required.
    - **Type**: `str`
    - **Example**: `-n HLA1_chr6`

# Prepare reference panel

First, we need to download the reference panel from the [UCSC server](https://cgl.gi.ucsc.edu/data/giraffe/construction/). The reference panel should be in the VCF format.


