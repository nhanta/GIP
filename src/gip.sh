#!/bin/bash

# Default values
num_cpus=0
num_gpus=0

# Function to display usage
usage() {
    echo "Usage: $0 -data_path <path> -data_name <name> -method <GIP|GAIN> -num_al <number> [-num_cpus <number>] [-num_gpus <number>]"
    exit 1
}

# Parse command line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -data_path) data_path="$2"; shift ;;
        -data_name) data_name="$2"; shift ;;
        -method) method="$2"; shift ;;
        -num_al) num_al="$2"; shift ;;
        -num_cpus) num_cpus="$2"; shift ;;
        -num_gpus) num_gpus="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; usage ;;
    esac
    shift
done

# Check if required parameters are provided
if [ -z "$data_path" ] || [ -z "$data_name" ] || [ -z "$method" ] || [ -z "$num_al" ]; then
    usage
fi

# Log the input parameters
echo "Running with the following parameters:"
echo "Data path: $data_path"
echo "Data name: $data_name"
echo "Method: $method"
echo "Number of alleles: $num_al"
echo "CPUs: $num_cpus"
echo "GPUs: $num_gpus"

# Create output directory
output=${data_path}/output/${method}
mkdir -p "${output}"
log=$(readlink -f "${output}/${data_name}".log)
echo "Log file: ${log}"

echo "Remove output files if existing:"| tee -a "${log}"
if [ -e "${output}"/"${data_name}".imputed."${method}".hap ]
then
  echo "- rm ${output}/${data_name}.imputed.${method}.hap" | tee -a "${log}"
  rm "${output}"/"${data_name}".imputed."${method}".hap
fi
if [ -e "${data_path}"/"${data_name}".hap ]
then
  echo "- rm ${data_path}/${data_name}.hap" | tee -a "${log}"
  rm "${data_path}"/"${data_name}".hap
fi
if [ -e "${data_path}"/"${data_name}".legend ]
then
  echo "- rm ${data_path}/${data_name}.legend" | tee -a "${log}"
  rm "${data_path}"/"${data_name}".legend
fi
if [ -e "${output}"/"${data_name}".imputed."${method}".vcf ]
then
  echo "- rm ${output}/${data_name}.imputed.${method}.vcf" | tee -a "${log}"
  rm "${output}"/"${data_name}".imputed."${method}".vcf
fi
if [ -e "${output}"/"${data_name}".imputed."${method}".dip ]
then
  echo "- rm ${output}/${data_name}.imputed.${method}.dip" | tee -a "${log}"
  rm "${output}"/"${data_name}".imputed."${method}".dip
fi

echo "Convert diplotype to haplotype: " | tee -a "${log}"
echo "- bcftools query -f '[%GT ]' ${data_path}/${data_name}.vcf | sed -E 's/[|\/]/ /g' | sed -E 's/\./?/g' | sed 's/ *$//' > ${data_path}/${data_name}.hap" | tee -a "${log}"
bcftools query -f '[%GT ]' "${data_path}"/"${data_name}".vcf | sed -E 's/[|\/]/ /g' | sed -E 's/\./?/g' | sed 's/ *$//' > "${data_path}"/"${data_name}".hap

echo "Create legend file: " | tee -a "${log}"
echo "- bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT' ${data_path}/${data_name}.vcf  > ${data_path}/${data_name}.legend" | tee -a "${log}"
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT' "${data_path}"/"${data_name}".vcf  > "${data_path}"/"${data_name}".legend

echo "Run imputation:" | tee -a "${log}"
echo """
python3 gip/main_new.py
--missing_data ${data_path}/${data_name}.hap
--method ${method}
--num_al ${num_al}
--output_data ${output}/${data_name}.imputed.${method}.hap
--batch_size 10
--hint_rate 0.9
--alpha 100
--iterations 10000
--num_cpus ${num_cpus}
--num_gpus ${num_gpus}
""" | tee -a "${log}"
SECONDS=0
python3 gip/main.py \
--missing_data "${data_path}"/"${data_name}".hap \
--method "${method}" \
--num_al "${num_al}" \
--output_data "${output}"/"${data_name}".imputed."${method}".hap \
--batch_size 10 --hint_rate 0.9 --alpha 100 \
--iterations 10000 \
--num_cpus "${num_cpus}" \
--num_gpus "${num_gpus}"
echo -e "- Sample: ${data_name}\n $(date -d@${SECONDS} -u +%H:%M:%S) seconds elapsed.\n" | tee -a "${log}"

echo "Export imputed data:" | tee -a "${log}"

echo "- Get header from VCF file: " | tee -a "${log}"
echo "  bcftools view -h ${data_path}/${data_name}.vcf | head -n -1 > ${output}/${data_name}.imputed.${method}.vcf" | tee -a "${log}"
bcftools view -h "${data_path}"/"${data_name}".vcf | head -n -1 > ${output}/${data_name}.imputed.${method}.vcf
echo -ne "#CHROM\tPOS\tID\tREF\tALT\t" >> ${output}/${data_name}.imputed.${method}.vcf
bcftools query -l ../data/HLA_3_alleles/HLA1_chr6.AFR.0.25.missing.vcf | tr '\n' '\t' >> ${output}/${data_name}.imputed.${method}.vcf

echo "- Convert imputed haplotype to diplotype: " | tee -a "${log}"
echo "  cat ${output}/${data_name}.imputed.${method}.hap | sed -E 's/([0-9?]*)\s([0-9?]*)\s/\1|\2\t/g' > ${output}/${data_name}.imputed.${method}.dip" | tee -a "${log}"
cat ${output}/${data_name}.imputed.${method}.hap | sed -E 's/([0-9?]*)\s([0-9?]*)\s/\1|\2\t/g' > ${output}/${data_name}.imputed.${method}.dip

echo "- Add legend and diplotype to file: " | tee -a "${log}"
echo "  paste -d'\t' ${data_path}/${data_name}.legend ${output}/${data_name}.imputed.${method}.dip >> ${output}/${data_name}.imputed.${method}.vcf" | tee -a "${log}"
paste -d'\t' ${data_path}/${data_name}.legend ${output}/${data_name}.imputed.${method}.dip >> ${output}/${data_name}.imputed.${method}.vcf