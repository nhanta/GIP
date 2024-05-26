import numpy as np
import pandas as pd
from utils import read_vcf, haploid_to_diploid, diploid_to_haploid
import argparse

def load_data (args):
    ori_data = args.ori_data
    missing_data = args.missing_data
    data_path = args.data_path

    # Load original data
    data = read_vcf(data_path)
    geno = data.iloc[:, 9::]
    col = geno.columns
    geno.rename(columns={col[-1]:col[-1].split("\n")[0]}, inplace=True)
    
    # Convert diploid to haploid
    hap = 

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
  load_data(args)