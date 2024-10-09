import pandas as pd
import numpy as np

def binary_sampler(p, rows, cols):
  '''Sample binary random variables.

  Args:
    - p: probability of 1
    - rows: the number of rows
    - cols: the number of columns

  Returns:
    - binary_random_matrix: generated binary random matrix.
  '''
  np.random.seed(7)
  unif_random_matrix = np.random.uniform(0., 1., size = [rows, cols])
  binary_random_matrix = 1*(unif_random_matrix < p)
  return binary_random_matrix

# Read a vcf file
def read_vcf(vcf_path):
  header=''
  with open(vcf_path, "rt") as ifile:
    for line in ifile:
      if line.startswith("##"):
        header += line
      elif line.startswith("#CHROM"):
        vcf_names = [x for x in line.split('\t')]
        break
  ifile.close()
  data = pd.read_csv(vcf_path, comment='#', sep="\s+", header=None, names=vcf_names)
  return header, data

def get_data(X, miss_rate):
  # Parameters
  no, dim = X.shape

  # Introduce missing data
  data_m = binary_sampler(1-miss_rate, no, dim)
  miss_data_x = X.copy()
  miss_data_x[data_m == 0] = ".|."
  return X, miss_data_x

def save_vcf(data, output_VCF, header = """##fileformat=VCFv4.1\n"""):
  with open(output_VCF, 'w') as vcf:
    vcf.write(header)
  data.to_csv(output_VCF, sep="\t", mode='a', index=False)