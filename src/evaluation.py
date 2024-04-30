import numpy as np
import pandas as pd
from metrics import eval_acc
from utils import read_vcf

ori_data = read_vcf("../data/HLA/HLA1_chr6.vcf").iloc[:, 9::]
imputed_data = read_vcf("../data/HLA/imputed_biallelic_geno.vcf").iloc[:, 9::]

mask = pd.read_csv("../data/HLA/mask.csv").set_index('Unnamed: 0')
evaluation = eval_acc(np.ravel(ori_data[mask==0]), np.ravel(imputed_data[mask==0]))
print(evaluation)
