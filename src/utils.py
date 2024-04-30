import numpy as np
import pandas as pd

# Read a vcf file
def read_vcf(vcf_path):
    with open(vcf_path, "rt") as ifile:
          for line in ifile:
            if line.startswith("#CHROM"):
                  vcf_names = [x for x in line.split('\t')]
                  break
    ifile.close()
    data = pd.read_csv(vcf_path, comment='#', delim_whitespace=True, header=None, names=vcf_names)
    return data
