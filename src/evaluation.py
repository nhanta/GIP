import pandas as pd
from metrics import eval_acc
import argparse

parser = argparse.ArgumentParser(
                    prog='GIP',
                    description='Data generation',
                    epilog='Help')

parser.add_argument('ori_data', metavar='n', type=str, help='original data')
parser.add_argument('missing_data', metavar='n', type=str, help='missing data')
parser.add_argument('imputed_data', metavar='r', type=str, help='imputed data')
args = parser.parse_args()

# Load original data and imputed data
ori_data = pd.read_csv(args.ori_data, sep = " ", header = None).to_numpy()
missing_data = pd.read_csv(args.missing_data, sep = " ", header = None).to_numpy()
#imputed_data = pd.read_csv(args.imputed_data, sep = "\t", header = None).to_numpy()
imputed_data = pd.read_csv(args.imputed_data, sep = " ", header = None).to_numpy()
print(ori_data)
print(missing_data)
print(imputed_data)
# Load mask data
data_m = missing_data == "?"

# Evaluate models
evaluation = eval_acc(ori_data[data_m], imputed_data[data_m])

