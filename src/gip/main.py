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

'''Main function for UCI letter and spam datasets.
'''

# Necessary packages
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import argparse
import numpy as np
import pandas as pd

from gip import imputation
from utils import save_hap, read_vcf, save_vcf, diploid_to_haploid, haploid_to_diploid

def main (args):
  '''Main function for UCI letter and spam datasets.
  
  Args:
    - data_name: letter or spam
    - miss_rate: probability of missing components
    - batch:size: batch size
    - hint_rate: hint rate
    - alpha: hyperparameter
    - iterations: iterations
    
  Returns:
    - imputed_data_x: imputed data
    - rmse: Root Mean Squared Error
  '''
  
  method = args.method
  missing_data = args.missing_data
  output_data = args.output_data
  num_al = args.num_al
  num_threads = args.num_threads
  
  gain_parameters = {'batch_size': args.batch_size,
                     'hint_rate': args.hint_rate,
                     'alpha': args.alpha,
                     'iterations': args.iterations}
  if num_al == 2:
    # Load missing data
    miss_data_x = pd.read_csv(missing_data, header = None, sep = " ").to_numpy()
    data_m = miss_data_x == "?"
  else:
    # Load original data
    data = read_vcf(missing_data[0:-3] + "vcf")
    geno = data.iloc[:, 9::]
    col = geno.columns
    geno.rename(columns={col[-1]:col[-1].split("\n")[0]}, inplace=True)
    
    # Convert diploid to haploid
    miss_data_x = diploid_to_haploid(geno).to_numpy()
    save_hap(miss_data_x, missing_data[0:-3] + "hap")
    data_m = miss_data_x == "."

  miss_data_x[data_m] = np.nan
  miss_data_x = miss_data_x.astype(float)
  # Impute missing data
  imputed_data_x = imputation(miss_data_x, gain_parameters, num_al).impute(method, num_threads)
  
  save_hap(imputed_data_x.astype(int).astype(str), output_data)

  if num_al != 2:
    ori_data = read_vcf(missing_data[0:-11] + "ori.vcf")
    ori_geno = ori_data.iloc[:, 9::]
    ori_col = geno.columns
    ori_geno.rename(columns={ori_col[-1]:ori_col[-1].split("\n")[0]}, inplace=True)
    
    # Convert diploid to haploid for original data
    ori_hap = diploid_to_haploid(ori_geno).to_numpy()
    save_hap(ori_hap, missing_data[0:-11] + "ori.hap")

    # Output the imputed diploid data
    im_data = pd.DataFrame(data = imputed_data_x.astype(int).astype(str), columns = list(range(2*len(ori_col))))
    save_hap(im_data, output_data[0:-3] + "imputed." + method + ".hap")
    imputed_data_value = haploid_to_diploid(im_data).to_numpy()
    imputed_data = pd.DataFrame(data = imputed_data_value, columns = ori_col)

    with open(missing_data[0:-3] + "header.txt") as f:
        # Read the contents of the file into a variable
        header = f.read()

    save_vcf(pd.concat([ori_data.iloc[:, 0:9], imputed_data], axis = 1), output_data[0:-3] + "imputed" + "." + method + ".vcf", header)


if __name__ == '__main__':  
  
  # Inputs for the main function
  parser = argparse.ArgumentParser()
  parser.add_argument(
      '--method',
      choices=['gip', 'gain'],
      default='gip',
      type=str)
  parser.add_argument(
      '--missing_data',
      help='missing data',
      default='gip',
      type=str)
  parser.add_argument(
      '--output_data',
      help='imputed data',
      default='gip',
      type=str)
  parser.add_argument(
      '--num_al',
      help='number of alleles',
      default=int,
      type=int)
  parser.add_argument(
      '--batch_size',
      help='the number of samples in mini-batch',
      default=128,
      type=int)
  parser.add_argument(
      '--hint_rate',
      help='hint probability',
      default=0.9,
      type=float)
  parser.add_argument(
      '--alpha',
      help='hyperparameter',
      default=100,
      type=float)
  parser.add_argument(
      '--iterations',
      help='number of training interations',
      default=10000,
      type=int)
  parser.add_argument(
      '--num_threads',
      help='number of training interations',
      default=20,
      type=int)
  
  args = parser.parse_args() 
  
  # Calls main function  
  main(args)