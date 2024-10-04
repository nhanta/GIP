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
output=${data_path}/output
mkdir -p "${output}"
# Remove output files if existing
[ -e "${output}"/"${data_name}".imputed."${method}".hap ] && rm "${output}"/"${data_name}".imputed."${method}".hap
[ -e "${data_path}"/"${data_name}".hap ] && rm "${data_path}"/"${data_name}".hap
[ -e "${data_path}"/"${data_name:0:-7}"ori.hap ] && rm "${data_path}"/"${data_name:0:-7}"ori.hap

gip() {
    log=$(readlink -f "${data_name}".log.txt)
    # Run imputation
    echo "python3 gip/main.py --missing_data "${data_path}"/"${data_name}".hap --method "${method}" --num_al "${num_al}" --output_data "${output}"/"${data_name}".imputed."${method}".hap --batch_size 10 --hint_rate 0.9 --alpha 100 --iterations 10000 --num_cpus "${num_cpus}" --num_gpus "${num_gpus}
    SECONDS=0
    python3 gip/main.py --missing_data "${data_path}"/"${data_name}".hap \
    --method "${method}" \
    --num_al "${num_al}" \
    --output_data "${output}"/"${data_name}".imputed."${method}".hap \
    --batch_size 10 --hint_rate 0.9 --alpha 100 \
    --iterations 10000 \
    --num_cpus "${num_cpus}" \
    --num_gpus "${num_gpus}"
    duration=$SECONDS
    echo -e "Sample: ${data_name}\n$((duration / 60)) minutes and $((duration % 60)) seconds elapsed.\n" | tee -a "${log}"
}

if [[ ${num_al} == 2 ]]; then
    # Diplotype to haplotype
    bcftools convert --haplegendsample "${data_path}"/"$data_name" "${data_path}"/"${data_name}".vcf
    gunzip -f "${data_path}"/"${data_name}".hap.gz
    gip
    # Export imputed data
    bgzip -f "${output}"/"${data_name}".imputed."${method}".hap
    cp "${data_path}"/"${data_name}".samples  "${output}"/"${data_name}".imputed."${method}".samples
    cp "${data_path}"/"${data_name}".legend.gz  "${output}"/"${data_name}".imputed."${method}".legend.gz

    bcftools convert --haplegendsample2vcf "${output}"/"${data_name}".imputed."${method}" > "${output}"/"${data_name}".imputed."${method}".vcf
else
    # Generate header
    grep '^##'  "${data_path}"/"$data_name".vcf > "${data_path}"/"$data_name".header.txt
    gip
fi

echo "Execution complete."