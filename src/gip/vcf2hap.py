import pandas as pd
from utils import read_vcf, save_hap, diploid_to_haploid
import argparse

parser = argparse.ArgumentParser(
                    prog='GIP',
                    description='Data generation',
                    epilog='Help')

parser.add_argument('ori_data', metavar='oi', type=str, help='original data')
parser.add_argument('missing_data', metavar='mi', type=str, help='missing data')
parser.add_argument('imputed_data', metavar='i', type=str, help='imputed data')
parser.add_argument('output_data', metavar='o', type=str, help='output data')
parser.add_argument('method', metavar='m', type=str, help='methods')
args = parser.parse_args()

# Load imputed data
im_data = read_vcf(args.imputed_data)

# Convert imputed data to haploid data
save_hap(diploid_to_haploid(im_data.iloc[:, 9::]).to_numpy(), args.output_data)

if args.method == "impute5":
    # Load original data
    ori = read_vcf(args.ori_data)
    ori.set_index("ID", inplace = True)

    missing = read_vcf(args.missing_data)
    missing.set_index("ID", inplace = True)

    # Get variants including in the original data
    inter = list(set(ori.index).intersection(set(im_data["ID"])))
    ori = ori.loc[inter, :]
    missing = missing.loc[inter, :]
    # Save the original haploid data for Impute5
    save_hap(diploid_to_haploid(ori.iloc[:, 8::]).to_numpy(), args.ori_data[0:-3] + "ip5.hap")
    # Save missing haploid data for Impute5
    save_hap(diploid_to_haploid(missing.iloc[:, 8::]).to_numpy(), args.missing_data[0:-3] + "ip5.hap")

    