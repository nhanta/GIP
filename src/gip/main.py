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

'''Main function for the 1000 Genomes Project (1KGP) phase 3.
'''

# Necessary packages
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import argparse
import numpy as np
import pandas as pd

from gip import GIP
from utils import save_hap

def main (args):
  '''Main function for the 1KGP phase 3.
  
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
  num_alleles = args.num_alleles

  gip_parameters = {'batch_size': args.batch_size,
                     'hint_rate': args.hint_rate,
                     'alpha': args.alpha,
                     'iterations': args.iterations,
                     'num_alleles': args.num_alleles
                     }
  system_parameters = {'num_cpus': args.num_cpus,
                       'num_gpus': args.num_gpus
                     }

  # Load missing data
  miss_data_x = pd.read_csv(missing_data, header = None, sep = " ").to_numpy()
  data_m = miss_data_x == "?"


  miss_data_x[data_m] = np.nan
  miss_data_x = miss_data_x.astype(float)
  # Impute missing data
  imputed_data_x = GIP(miss_data_x, gip_parameters, system_parameters).impute(method)

  save_hap(imputed_data_x.astype(int).astype(str), output_data)

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
      '--num_alleles',
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
      '--num_cpus',
      help='number of cpus',
      default=2,
      type=int)
  parser.add_argument(
      '--num_gpus',
      help='number of gpus',
      default=0,
      type=int)

  args = parser.parse_args() 
  
  # Calls main function  
  main(args)