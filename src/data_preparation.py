import pandas as pd
from utils import get_data, read_vcf, save_vcf
import argparse

parser = argparse.ArgumentParser(
                    prog='GIP',
                    description='Data generation',
                    epilog='Help')

parser.add_argument('data_name', metavar='n', type=str, help='Prefix of data name')
parser.add_argument('rate', metavar='r', type=str, help='Missing rate')
args = parser.parse_args()

rate = float(args.rate)
vcf_path = "../data/HLA/" + args.data_name + ".vcf"
output_VCF = "../data/HLA/" + args.data_name + "." + str(rate) + ".missing.vcf"
output_mask = "../data/HLA/" + args.data_name + "." + str(rate) + ".mask.csv"


if __name__ == '__main__':

    print("")
    print("")
    print("|============================================================================|")
    print("|                                                                            |")
    print("|         -----            DATA PREPARATION FOR GIP          -----           |")
    print("|                                                                            |")
    print("|============================================================================|")
    print("")
    print("")
    print("********************************* INPUT DATA *********************************")
    print("")
    print("Import data may take several minutes, please wait...")
    print("")

    # Load original data
    data = read_vcf(vcf_path)
    geno = data.iloc[:, 9::]
    geno, geno_miss, m = get_data(geno, rate)
    gn = pd.concat([data.iloc[:, 0:9], geno_miss], axis = 1)
    gn.rename(columns={'NA21144\n':'NA21144'}, inplace=True)
    
    print("********************************** SAVING **********************************")
    save_vcf(gn, output_VCF)
    pd.DataFrame(columns = list(gn.columns)[9::], data = m).to_csv(output_mask)
    print("")
    print("********************************* FINISHED *********************************")
    print("")
    