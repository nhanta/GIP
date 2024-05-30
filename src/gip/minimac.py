import pandas as pd
from utils import read_vcf, save_hap, diploid_to_haploid
import argparse

parser = argparse.ArgumentParser(
                    prog='GIP',
                    description='Data generation',
                    epilog='Help')

#parser.add_argument('ori_data', metavar='oi', type=str, help='original data')
parser.add_argument('imputed_data', metavar='i', type=str, help='imputed data')
parser.add_argument('output_data', metavar='o', type=str, help='output data')
args = parser.parse_args()

'''
# Load original data
ori = read_vcf(args.ori_data)
ori.set_index("ID", inplace = True)

# Load imputed data
im_data = read_vcf(args.imputed_data)
im_data.set_index("ID", inplace = True)

# Get variants including in the original data
inter = set(ori.index).intersection(set(im_data))
imp_data = im_data.iloc[list(inter), :]
imp_data.reset_index(inplace=True)
print(imp_data)
'''
# Load imputed data
im_data = read_vcf(args.imputed_data)

# Convert imputed data to haploid data
save_hap(diploid_to_haploid(im_data.iloc[:, 9::]).to_numpy(), args.output_data)
