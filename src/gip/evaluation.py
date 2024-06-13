import numpy as np
import pandas as pd
from metrics import eval_acc
import argparse

parser = argparse.ArgumentParser(
                    prog='GIP',
                    description='Data generation',
                    epilog='Help')

parser.add_argument('ori_data', metavar='oi', type=str, help='original data')
parser.add_argument('missing_data', metavar='mi', type=str, help='missing data')
parser.add_argument('imputed_data', metavar='ii', type=str, help='imputed data')
parser.add_argument('num_al', metavar='n', type=int, help='Number of alleles')
parser.add_argument('method', metavar='m', type=str, help='Method')
args = parser.parse_args()

# Load original data and imputed data
ori_data = pd.read_csv(args.ori_data, sep = "\s+", header = None).to_numpy().astype(int)
missing_data = pd.read_csv(args.missing_data, sep = "\s+", header = None, low_memory=False).to_numpy()
imputed_data = pd.read_csv(args.imputed_data, sep = "\s+", header = None).to_numpy().astype(int)

# Load mask data

if np.sum(missing_data == "?") != 0:
    data_m = missing_data == "?"
else:
    data_m = missing_data == "."

# Evaluate models
# evaluation = eval_acc(ori_data[data_m], imputed_data[data_m])
evaluation = eval_acc(ori_data.flatten(), imputed_data.flatten())

name = args.ori_data.split("/")[-1][0:-7] + args.method + "." + str(args.num_al)
with open("results." + str(args.num_al)  + "_alleles.txt", 'a+') as f:
     f.write(name + ":" + str(evaluation) + '\n')
     f.write("")
     f.close()

