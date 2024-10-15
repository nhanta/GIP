import numpy as np
import pandas as pd
import argparse

if __name__ == 'main':
  parser = argparse.ArgumentParser(
    prog='GIP',
    description='Data evaluation',
    epilog='Help')
  parser.add_argument('-o','--ori_data', metavar='oi', type=str, help='Original data (hap)')
  parser.add_argument('-m','--missing_data', metavar='mi', type=str, help='Missing data (hap)')
  parser.add_argument('-i','--imputed_data', metavar='ii', type=str, help='Imputed data (hap)')
  parser.add_argument('-p','--path', metavar='n', type=int, help='Result path')
  parser.add_argument('-m','--method', metavar='m', type=str, help='Method')
  args = parser.parse_args()
  # Load original data and imputed data
  ori_data = pd.read_csv(args.ori_data, sep = "\s+", header = None).to_numpy().astype(int)
  imputed_data = pd.read_csv(args.imputed_data, sep = "\s+", header = None).to_numpy().astype(int)
  ori_flat = ori_data.flatten()
  imp_flat = imputed_data.flatten()
  acc = sum(ori_flat == imp_flat) / len(ori_flat)
  print(f"Accuracy: {acc}")
  with open(args.path, 'a+') as f:
    f.write(f"""
{args.method}: {acc}
\torigin: {args.ori_data}
\tmissing: {args.missing_data}
\timputed: {args.imputed_data}
""")
    f.close()
