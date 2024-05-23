#!/bin/bash
data_path=$1 
data_name=$2 # Prefix name of data
method=$3 
num_al=$4
num_threads=$5

# Create output directory
output=${data_path}/output
mkdir ${output}
parser.add_argument('ori_data', metavar='n', type=str, help='original data')
parser.add_argument('missing_data', metavar='n', type=str, help='missing data')
parser.add_argument('imputed_data', metavar='r', type=str, help='imputed data')
gunzip ${output}/${data_name}.imputed.hap.gz