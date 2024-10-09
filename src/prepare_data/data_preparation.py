import pandas as pd
from pathlib import Path
from utils import get_data, read_vcf, save_vcf
import argparse
if __name__ == '__main__':
  # sp_code = ['all', 'AFR', 'AMR', 'EAS', 'EUR', 'SAS']
  parser = argparse.ArgumentParser(
    prog='GIP',
    description='Data generation',
    epilog='Help')
  parser.add_argument('-p', '--data_path', metavar='p', type=str,
                      help='Data path containing vcf and metadata')
  parser.add_argument('-s', '--superpopulation_code', metavar='s', default=['all'],
                      help='Super population code will be selected', type=str, nargs='+')
  parser.add_argument('-r', '--rate', metavar='r', type=str, default=['0.25'],
                      help='Rate of missing data', nargs = '+')
  parser.add_argument('-m', '--meta', metavar='m', type=str, help='Metadata file')
  parser.add_argument('-n', '--data_name', metavar='n', type=str, help='Prefix of data name')
  args = parser.parse_args()

  print("""
|============================================================================|
|                                                                            |
|         -----            DATA PREPARATION FOR GIP          -----           |
|                                                                            |
|============================================================================|

********************************* INPUT DATA *********************************

Import data may take several minutes, please wait...

""")

  vcf_path = Path(args.data_path, args.data_name + ".vcf")
  meta_path = Path(args.data_path, args.meta)

  # Load data
  header, data = read_vcf(vcf_path)
  data.columns = [c.replace("\n", "_") for c in data.columns]
  geno = data.iloc[:, 9::]
  pop = geno.columns
  print("Number of all population: ", len(list(pop)))
  print("Number of positions: ", geno.shape[0])

  # Load meta_data
  meta_data = pd.read_csv(meta_path, sep="\t")
  # split sample by meta data
  geno_data = {'all': geno}
  for code in args.superpopulation_code:
    if code == 'all': continue
    sample = meta_data["Sample name"][meta_data["Superpopulation code"] == code]
    inter = set(sample).intersection(set(pop))
    print("Number of ", code, ": ", len(list(inter)))
    geno_data[code] = geno.loc[:, list(inter)]

  print("""
********************************** SAVING **********************************
""")
  for p in args.superpopulation_code:
    print("Population: ", p)
    for rate in args.rate:
      print("\tRate: ", rate)
      output_miss_path = Path(args.data_path, args.data_name + "." + p + "." + str(rate) + ".missing.vcf")
      output_ori_data = Path(args.data_path, args.data_name + "." + p + "." + str(rate) + ".ori.vcf")
      gn, geno_miss = get_data(geno_data[p], float(rate))
      miss_data = pd.concat([data.iloc[:, 0:9], geno_miss], axis=1)
      ori_data = pd.concat([data.iloc[:, 0:9], gn], axis=1)
      save_vcf(miss_data, output_miss_path, header)
      save_vcf(ori_data, output_ori_data, header)
      print("\t\tSaved missing data: ", output_miss_path)
      print("\t\tSaved original data: ", output_ori_data)
print("""
********************************* FINISHED *********************************
""")