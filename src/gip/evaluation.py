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
args = parser.parse_args()

# Load original data and imputed data
ori_data = pd.read_csv(args.ori_data, sep = "\t", header = None).to_numpy()
missing_data = pd.read_csv(args.missing_data, sep = "\t", header = None).to_numpy()
imputed_data = pd.read_csv(args.imputed_data, sep = "\t", header = None).to_numpy()

# Load mask data
if args.num_al == 2:
    data_m = missing_data == "?"
else:
    data_m = missing_data == "."

# Evaluate models
evaluation = eval_acc(ori_data[data_m], imputed_data[data_m])

