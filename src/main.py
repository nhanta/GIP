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
from utils import rmse_loss
from metrics import eval_acc

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
  ori_data = args.ori_data
  output_data = args.output_data
  mask_data = args.mask_data
  num_threads = args.num_threads
  
  gain_parameters = {'batch_size': args.batch_size,
                     'hint_rate': args.hint_rate,
                     'alpha': args.alpha,
                     'iterations': args.iterations}
  
  # Load original data and imputed data
  ori_data_x = pd.read_csv(ori_data, sep = "\s+", header = None)
  ori_data_x = ori_data_x.iloc[:, 5::].to_numpy()
  ori_data_x = ori_data_x.astype(float)
  
  # Load mask data
  data_m  = pd.read_csv(mask_data).set_index('Unnamed: 0')
  data_m = pd.concat([data_m, data_m], axis = 1) # Convert diploid mask to haploid mask
  data_m = data_m.to_numpy()
  
  # Missing data
  miss_data_x = ori_data_x.copy()
  miss_data_x[data_m == 0] = np.nan
  
  # Impute missing data
  imputed_data_x = imputation(miss_data_x, gain_parameters).impute(method, num_threads)

  # Report the RMSE performance
  rmse = rmse_loss (ori_data_x, imputed_data_x, data_m)
  acc = eval_acc(np.ravel(ori_data_x[data_m==0]), np.ravel(imputed_data_x[data_m==0]))
  print()
  print('RMSE Performance: ' + str(np.round(rmse, 4)))
  #evaluation.to_csv("results/gip_evaluated_" + data_name)
  imputed_data_x.to_csv(output_data, header = False, index = False)
  return imputed_data_x, acc, rmse

if __name__ == '__main__':  
  
  # Inputs for the main function
  parser = argparse.ArgumentParser()
  parser.add_argument(
      '--method',
      choices=['gip', 'gain'],
      default='gip',
      type=str)
  parser.add_argument(
      '--ori_data',
      help='original data',
      default='gip',
      type=str)
  parser.add_argument(
      '--output_data',
      help='imputed data',
      default='gip',
      type=str)
  parser.add_argument(
      '--mask_data',
      help='mask data',
      default='spam',
      type=str)
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
  imputed_data, acc, rmse = main(args)