import numpy as np
import pandas as pd
from metrics import eval
from utils import read_vcf

# Load original data and imputed data
ori_data = pd.read_csv("../data/HLA/HLA1_chr6.haps", sep = "\s+", header = None).iloc[:, 5::]
imputed_data = pd.read_csv("../data/HLA/imputed_biallelic_geno.haps", sep = "\s+", header = None).iloc[:, 5::]
# Load mask data
mask = pd.read_csv("../data/HLA/HLA1_chr6.0.2.mask.csv").set_index('Unnamed: 0')
mask = pd.concat([mask, mask], axis = 1) # Convert diploid mask to haploid mask
print(mask.shape)
ori_data = ori_data.to_numpy()
imputed_data = imputed_data.to_numpy()
mask = mask.to_numpy()

# Evaluate models
evaluation = eval(np.ravel(ori_data[mask==0]), np.ravel(imputed_data[mask==0]))
