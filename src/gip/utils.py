# coding=utf-8
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

'''Utility functions for GIP.

(1) normalization: MinMax Normalizer
(2) renormalization: Recover the data from normalzied data
(3) rounding: Handlecategorical variables after imputation
(4) rmse_loss: Evaluate imputed data in terms of RMSE
(5) binary_sampler: sample binary random variables
(6) uniform_sampler: sample uniform random variables
(7) sample_batch_index: sample random batch index
'''
 
# Necessary packages
import numpy as np
import pandas as pd

def normalization (data, parameters=None):
  '''Normalize data in [0, 1] range.
  
  Args:
    - data: original data
  
  Returns:
    - norm_data: normalized data
    - norm_parameters: min_val, max_val for each feature for renormalization
  '''

  # Parameters
  _, dim = data.shape
  norm_data = data.copy()

  if parameters is None:
  
    # MixMax normalization
    min_val = np.zeros(dim)
    max_val = np.zeros(dim)
    
    # For each dimension
    for i in range(dim):
      min_val[i] = np.nanmin(norm_data[:,i])
      norm_data[:,i] = norm_data[:,i] - np.nanmin(norm_data[:,i])
      max_val[i] = np.nanmax(norm_data[:,i])
      norm_data[:,i] = norm_data[:,i] / (np.nanmax(norm_data[:,i]) + 1e-6)   
      
    # Return norm_parameters for renormalization
    norm_parameters = {'min_val': min_val,
                       'max_val': max_val}

  else:
    min_val = parameters['min_val']
    max_val = parameters['max_val']
    
    # For each dimension
    for i in range(dim):
      norm_data[:,i] = norm_data[:,i] - min_val
      norm_data[:,i] = norm_data[:,i] / (max_val + 1e-6)  
      
    norm_parameters = parameters    
      
  return norm_data, norm_parameters


def renormalization (norm_data, norm_parameters):
  '''Renormalize data from [0, 1] range to the original range.
  
  Args:
    - norm_data: normalized data
    - norm_parameters: min_val, max_val for each feature for renormalization
  
  \Returns:
    - renorm_data: renormalized original data
  '''
  
  min_val = norm_parameters['min_val']
  max_val = norm_parameters['max_val']

  _, dim = norm_data.shape
  renorm_data = norm_data.copy()
  try:
    for i in range(dim):
      renorm_data[:,i] = renorm_data[:,i] * (max_val[i] + 1e-6)   
      renorm_data[:,i] = renorm_data[:,i] + min_val[i]
  except:
     for i in range(dim):
      renorm_data[:,i] = renorm_data[:,i] * (max_val + 1e-6)   
      renorm_data[:,i] = renorm_data[:,i] + min_val
    
  return renorm_data


def rounding (imputed_data, data_x):
  '''Round imputed data for categorical variables.
  
  Args:
    - imputed_data: imputed data
    - data_x: original data with missing values
    
  Returns:
    - rounded_data: rounded imputed data
  '''
  
  _, dim = data_x.shape
  rounded_data = imputed_data.copy()
  
  for i in range(dim):
    temp = data_x[~np.isnan(data_x[:, i]), i]
    # Only for the categorical variable
    if len(np.unique(temp)) < 20:
      rounded_data[:, i] = np.round(rounded_data[:, i])
      
  return rounded_data


def rmse_loss (ori_data, imputed_data, data_m):
  '''Compute RMSE loss between ori_data and imputed_data
  
  Args:
    - ori_data: original data without missing values
    - imputed_data: imputed data
    - data_m: indicator matrix for missingness
    
  Returns:
    - rmse: Root Mean Squared Error
  '''
  
  ori_data, norm_parameters = normalization(ori_data)
  imputed_data, _ = normalization(imputed_data, norm_parameters)
    
  # Only for missing values
  nominator = np.sum(((1-data_m) * ori_data - (1-data_m) * imputed_data)**2)
  denominator = np.sum(1-data_m)
  
  rmse = np.sqrt(nominator/float(denominator))
  
  return rmse

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


def uniform_sampler(low, high, rows, cols):
  '''Sample uniform random variables.
  
  Args:
    - low: low limit
    - high: high limit
    - rows: the number of rows
    - cols: the number of columns
    
  Returns:
    - uniform_random_matrix: generated uniform random matrix.
  '''
  np.random.seed(7)
  return np.random.uniform(low, high, size = [rows, cols])       


def sample_batch_index(total, batch_size):
  '''Sample index of the mini-batch.
  
  Args:
    - total: total number of samples
    - batch_size: batch size
    
  Returns:
    - batch_idx: batch index
  '''
  np.random.seed(7)
  total_idx = np.random.permutation(total)
  batch_idx = total_idx[:batch_size]
  return batch_idx
  
def prob(num_al, N):
  '''
  Calculate probability(h_{k+1,j}|h_1,...,h_k)

  Args:
    - X: a matrix with missing genotypes
    
  Returns:
    - p: a probability
  '''

  k = N - 1

  # Calculate theta
  # Number of alleles -1 (num_al - 1) is number of mutation in a locus
  theta = 0
  for m in range(1, k + 1):
    theta += 1/m
  theta = (num_al - 1)/theta
  
  # Probability 
  A = k/(k + theta) 
  # Probability
  B = (1/num_al) * (theta/(k + theta))
    
  return A, B
    

# Read a vcf file
def read_vcf(vcf_path):
    with open(vcf_path, "rt") as ifile:
          for line in ifile:
            if line.startswith("#CHROM"):
                  vcf_names = [x for x in line.split('\t')]
                  break
    ifile.close()
    data = pd.read_csv(vcf_path, comment='#', sep="\s+", header=None, names=vcf_names)
    return data

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

def save_hap(data, output_path):
   cols = list(range(data.shape[1]))
   df = pd.DataFrame(data = data, columns = cols)
   df.to_csv(output_path, header = False, index = False, sep ="\t", mode='a')

# Convert diploid to haploid
def diploid_to_haploid(data):    
# Define the columns to split and the separator
    columns_to_split = data.columns
    df = pd.DataFrame()
    value = []
    # Split the specified columns
    for col in columns_to_split:
        try:
            # Split the column into multiple columns
            split_cols = data[col].str.split('|', expand=True)
        except:
            split_cols = data[col].str.split('/', expand=True)

        # Rename the split columns to avoid name clashes
        split_cols.columns = [f"{col}_{i+1}" for i in range(split_cols.shape[1])]
        
        # Concatenate the split columns back to the original DataFrame
        df = pd.concat([df, split_cols], axis=1)
        
    return(df)

# Convert haploid to diploid with data frame
def haploid_to_diploid(df):    
    # Function to merge two columns with a separator
    def merge_columns(df, col1, col2):
        return df[col1].astype(str) + "|" + df[col2].astype(str)

    # List to hold the merged columns
    merged_columns = []
    l = df.shape[1] - 1
    
    # Loop through pairs of adjacent columns
    for i in range(0, l, 2):
        merged_col = merge_columns(df, df.columns[i], df.columns[i+1])
        merged_columns.append(merged_col)

    # Create a new DataFrame with the merged columns
    merged_df = pd.DataFrame(merged_columns).T

    return(merged_df)