import pandas as pd
from utils import get_data, read_vcf, save_vcf
import argparse

parser = argparse.ArgumentParser(
                    prog='GIP',
                    description='Data generation',
                    epilog='Help')

parser.add_argument('data_name', metavar='n', type=str, help='Prefix of data name')
args = parser.parse_args()
vcf_path = "../../data/HLA/" + args.data_name + ".vcf"

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
    data.rename(columns={'NA21144\n':'NA21144'}, inplace=True)
    geno = data.iloc[:, 9::]
    print("Number of all population: ", len(list(geno.columns)))

    # Load meta_data
    meta_data = pd.read_csv("../../data/HLA/igsr-1000 genomes on grch38.tsv", sep = "\t")

    # Select sub populations
    AFR = meta_data["Sample name"][meta_data["Superpopulation code"] == 'AFR']
    AMR = meta_data["Sample name"][meta_data["Superpopulation code"] == 'AMR']
    EAS = meta_data["Sample name"][meta_data["Superpopulation code"] == 'EAS']
    EUR = meta_data["Sample name"][meta_data["Superpopulation code"] == 'EUR']
    SAS = meta_data["Sample name"][meta_data["Superpopulation code"] == 'SAS']
    pop = geno.columns

    
    AFR_inter = set(AFR).intersection(set(pop)) 
    geno_afr = geno.loc[:, list(AFR_inter)]
    print ("Number of AFR: ", len(list(AFR_inter)))

    AMR_inter = set(AMR).intersection(set(pop)) 
    geno_amr = geno.loc[:, list(AMR_inter)]
    print ("Number of AMR: ", len(list(AMR_inter)))

    EAS_inter = set(EAS).intersection(set(pop)) 
    geno_eas = geno.loc[:, list(EAS_inter)]
    print ("Number of EAS: ", len(list(EAS_inter)))

    EUR_inter = set(EUR).intersection(set(pop)) 
    geno_eur = geno.loc[:, list(EUR_inter)]
    print ("Number of EUR: ", len(list(EUR_inter)))

    SAS_inter = set(SAS).intersection(set(pop)) 
    geno_sas = geno.loc[:, list(SAS_inter)]
    print ("Number of SAS: ", len(list(SAS_inter)))

    print("********************************** SAVING **********************************")

    for rate in [0.02, 0.05, 0.1, 0.2, 0.25]:
        for p in ['all', 'AFR', 'AMR', 'EAS', 'EUR', 'SAS']:
            output_miss_path = "../../data/HLA/" + args.data_name + "." + p \
                + "." + str(rate) + ".missing.vcf"
            output_ori_data = "../../data/HLA/" + args.data_name + "." + p \
                + "." + str(rate) + ".ori.vcf"
            if p == 'all':
                gn, geno_miss = get_data(geno, rate)
                miss_data = pd.concat([data.iloc[:, 0:9], geno_miss], axis = 1)
                ori_data = pd.concat([data.iloc[:, 0:9], gn], axis = 1)
                save_vcf(miss_data, output_miss_path)
                save_vcf(ori_data, output_ori_data)

            elif p == 'AFR':
                gn, geno_miss = get_data(geno_afr, rate)
                miss_data = pd.concat([data.iloc[:, 0:9], geno_miss], axis = 1)
                ori_data = pd.concat([data.iloc[:, 0:9], gn], axis = 1)
                save_vcf(miss_data, output_miss_path)
                save_vcf(ori_data, output_ori_data)

            elif p == 'AMR':
                gn, geno_miss = get_data(geno_amr, rate)
                miss_data = pd.concat([data.iloc[:, 0:9], geno_miss], axis = 1)
                ori_data = pd.concat([data.iloc[:, 0:9], gn], axis = 1)
                save_vcf(miss_data, output_miss_path)
                save_vcf(ori_data, output_ori_data)

            elif p == 'EAS':
                gn, geno_miss = get_data(geno_eas, rate)
                miss_data = pd.concat([data.iloc[:, 0:9], geno_miss], axis = 1)
                ori_data = pd.concat([data.iloc[:, 0:9], gn], axis = 1)
                save_vcf(miss_data, output_miss_path)
                save_vcf(ori_data, output_ori_data)
            
            elif p == 'EUR':
                gn, geno_miss = get_data(geno_eur, rate)
                miss_data = pd.concat([data.iloc[:, 0:9], geno_miss], axis = 1)
                ori_data = pd.concat([data.iloc[:, 0:9], gn], axis = 1)
                save_vcf(miss_data, output_miss_path)
                save_vcf(ori_data, output_ori_data)
            
            elif p == 'SAS':
                gn, geno_miss = get_data(geno_sas, rate)
                miss_data = pd.concat([data.iloc[:, 0:9], geno_miss], axis = 1)
                ori_data = pd.concat([data.iloc[:, 0:9], gn], axis = 1)
                save_vcf(miss_data, output_miss_path)
                save_vcf(ori_data, output_ori_data)

    print("")
    print("********************************* FINISHED *********************************")
    print("")
    