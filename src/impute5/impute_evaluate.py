import pandas as pd
import os
import glob
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=SyntaxWarning)
warnings.simplefilter(action='ignore', category=DeprecationWarning)

# config
missing_postfix = ".missing.vcf"
original_postfix = ".ori.vcf"
imputed_postfix = ".missing.vcf.final.vcf"

# Load data
os.chdir("impute5_test")
sampleList = glob.glob('*' + original_postfix)
samples = [sample.replace(original_postfix, '') for sample in sampleList]


# Read a vcf file
def read_vcf(vcf_path):
    with open(vcf_path, "rt") as ifile:
        for line in ifile:
            if line.startswith("#CHROM"):
                vcf_names = [x for x in line.split('\t')]
                break
    ifile.close()
    data = pd.read_csv(vcf_path, comment='#', sep="\s+", header=None, names=vcf_names)
    data['POS'] = data[['#CHROM', 'POS']].apply(lambda x: '{}.{}'.format(x[0], x[1]).replace("chr", ""), axis=1)
    data = data.drop(['#CHROM', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'], axis=1)
    return data


# validate
total_position_total = 0
missing_count_total = 0  # missing allele
wrong_with_missing_total = 0  # check if non-hidden allele change from missing file
wrong_with_original_total = 0  # check if hidden allele different from original file
wrong_position_total = 0  # check if it has position not in origin

for sample in samples:
    missing = read_vcf(sample + missing_postfix)
    imputed = read_vcf(sample + imputed_postfix)
    original = read_vcf(sample + original_postfix)
    true_position = list(original['POS'])

    wrong_with_missing = 0  # check if non-hidden allele change from missing file
    wrong_with_original = 0  # check if hidden allele different from original file
    wrong_position = 0  # check if it has position not in origin
    missing_count = 0  # check if it has position not in origin

    for i in range(imputed.shape[0]):
        if imputed.iloc[i]['POS'] not in true_position:
            wrong_position += 1
            continue
        position_in_original = true_position.index(imputed.iloc[i]['POS'])
        for j in range(1, imputed.shape[1]):
            if missing.iloc[position_in_original, j] == ".|." or missing.iloc[position_in_original, j] == "./.":
                missing_count += 2
                if original.iloc[position_in_original, j] != imputed.iloc[i, j]:
                    origin_diplo = original.iloc[position_in_original, j]
                    imputed_diplo = imputed.iloc[i, j]
                    if imputed_diplo[0] != origin_diplo[0]:
                        wrong_with_original += 1
                    if imputed_diplo[2] != origin_diplo[2]:
                        wrong_with_original += 1
            else:
                if missing.iloc[position_in_original, j] != imputed.iloc[i, j]:
                    missing_diplo = missing.iloc[position_in_original, j]
                    imputed_diplo = imputed.iloc[i, j]
                    if imputed_diplo[0] != missing_diplo[0]:
                        wrong_with_missing += 1
                    if imputed_diplo[2] != missing_diplo[2]:
                        wrong_with_missing += 1
    total_position = original.shape[1] * (original.shape[0] - 1) * 2
    print()
    warnings.warn("Done "+sample,Warning)
    print("Sample: {}".format(sample))
    print("Total allele: {}".format(total_position))
    print("Missing allele: {}".format(missing_count))
    print(
        "Change non-missing position in missing vcf: {:.2f}% ({})".format(
            wrong_with_missing / total_position * 100,
            wrong_with_missing))
    print(
        "Different from original with missing position: {:.2f}% ({})".format(
            wrong_with_original / total_position * 100,
            wrong_with_original))
    print("Create new position from original: {}".format(wrong_position))
    accuracy = (1 - (wrong_with_missing + wrong_with_original) / total_position) * 100
    print("Accuracy Total: {:.2f}%".format(accuracy))
    accuracy_on_missing = (1 - wrong_with_original / missing_count) * 100
    print("Accuracy On Missing: {:.2f}%".format(accuracy_on_missing))

    # add to total
    total_position_total += total_position
    missing_count_total += missing_count
    wrong_with_missing_total += wrong_with_missing
    wrong_with_original_total += wrong_with_original
    wrong_position_total += wrong_position

# print total info
print()
print("Total: {} samples".format(len(samples)))
print("Total position: {}".format(total_position_total))
print("Missing allele: {}".format(missing_count_total))
print("Change non-missing position in missing vcf: {:.2f}% ({})".format(
    wrong_with_missing_total / total_position_total * 100,
    wrong_with_missing_total))
print(
    "Different from original with missing position: {:.2f}% ({})".format(
        wrong_with_original_total / total_position_total * 100,
        wrong_with_original_total))
print("Create new position from original: {}".format(wrong_position_total))
accuracy_total = (1 - (wrong_with_missing_total + wrong_with_original_total) / total_position_total) * 100
print("Accuracy Total: {:.2f}%".format(accuracy_total))
accuracy_on_missing_total = (1 - wrong_with_original_total / missing_count_total) * 100
print("Accuracy On Missing: {:.2f}%".format(accuracy_on_missing_total))
