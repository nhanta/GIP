import pandas as pd
from utils import get_data, read_vcf, save_vcf
import argparse

parser = argparse.ArgumentParser(
                    prog='GIP',
                    description='Data generation',
                    epilog='Help')

parser.add_argument('data_path', metavar='p', type=str, help='Data path')
parser.add_argument('data_name', metavar='n', type=str, help='Prefix of data name')
args = parser.parse_args()

vcf_path = args.data_path + "/" + args.data_name + ".vcf"

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

    # Load header
    with open(args.data_path + "/" + "Chr1_DS.txt") as f:
        # Read the contents of the file into a variable
        header = f.read()

    # Load original data
    data = read_vcf(vcf_path)
    data.rename(columns={'NA21144\n':'NA21144'}, inplace=True)
    geno = data.iloc[:, 9::]
    print("Number of all population: ", len(list(geno.columns)))

    print("********************************** SAVING **********************************")

    for rate in [0.1, 0.2, 0.3]:
        for p in ['all']:
            output_miss_path = args.data_path + "/" + args.data_name + "." + p \
                + "." + str(rate) + ".missing.vcf"
            output_ori_data = args.data_path + "/" + args.data_name + "." + p \
                + "." + str(rate) + ".ori.vcf"
            if p == 'all':
                gn, geno_miss = get_data(geno, rate)
                miss_data = pd.concat([data.iloc[:, 0:9], geno_miss], axis = 1)
                ori_data = pd.concat([data.iloc[:, 0:9], gn], axis = 1)
                save_vcf(miss_data, output_miss_path, header)
                save_vcf(ori_data, output_ori_data, header)

    print("")
    print("********************************* FINISHED *********************************")
    print("")
    