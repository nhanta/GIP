beagle=$1 # java -Xmx96g -jar ../bin/beagle5.jar
beagle=${beagle:-"beagle"}
input_dir=$2 # ../data/HLA_3_alleles
file_name=$3 # HLA1_chr6.all.0.25.missing.vcf
output_dir=$4 # ../data/HLA_3_alleles/
ref_panel=$5 # ../data/other/ref_panel_multi_allele/1000GP_AF_chr6.vcf.gz
genetic_map=$6 # ../data/other/genetic_map/beagle_chr6_b38.map
num_threads=$7
num_threads=${num_threads:-$(nproc)}

echo "Running with the following parameters:"
echo "- Input directory: $input_dir"
echo "- File name: $file_name"
echo "- Output directory: $output_dir"
echo "- Reference panel: $ref_panel"
echo "- Genetic map: $genetic_map"
echo "- Number of threads: $num_threads"

mkdir -p ${output_dir}/output/beagle
output=${output_dir}/output/beagle
log=$(readlink -f "${output}/${file_name}".log)
echo "Log file: ${log}"

(
echo "Remove output files if existing:"
if [ -e "${output}"/"${file_name}".phased.vcf.gz ]
then
  echo "- rm ${output}/${file_name}.phased.vcf.gz"
  rm "${output}"/"${file_name}".phased.vcf.gz
fi
if [ -e "${output}"/"${file_name}".phased.log ]
then
  echo "- rm ${output}/${file_name}.phased.log"
  rm "${output}"/"${file_name}".phased.log
fi
if [ -e "${output}"/"${file_name}".imputed.vcf.gz ]
then
  echo "- rm ${output}/${file_name}.imputed.vcf.gz"
  rm "${output}"/"${file_name}".imputed.vcf.gz
fi
if [ -e "${output}"/"${file_name}".imputed.log ]
then
  echo "- rm ${output}/${file_name}.imputed.log"
  rm "${output}"/"${file_name}".imputed.log
fi
if [ -e "${output}"/"${file_name}".imputed.hap ]
then
  echo "- rm ${output}/${file_name}.imputed.hap"
  rm "${output}"/"${file_name}".imputed.hap
fi
echo "Phasing and imputing by Beagle5.4:"
echo "- $beagle gt=${input_dir}/${file_name} ref=${ref_panel} map=${genetic_map} out=${output}/${file_name}.phased impute=false nthreads=${num_threads}"
SECONDS=0
# Phasing by Beagle
$beagle \
  gt=${input_dir}/${file_name} \
  ref=${ref_panel} \
  map=${genetic_map} \
  out=${output}/${file_name}.phased \
  impute=false \
  gp=true \
  nthreads=${num_threads}
echo "- $beagle gt=${output}/${file_name}.phased.vcf.gz ref=${ref_panel} map=${genetic_map} out=${output}/${file_name}.imputed gp=true ap=true nthreads=${num_threads}"
# Impute by beagle
$beagle \
  gt=${output}/${file_name}.phased.vcf.gz \
  ref=${ref_panel} \
  map=${genetic_map} \
  out=${output}/${file_name}.imputed \
  gp=true \
  ap=true \
  nthreads=${num_threads}
echo -e "- Sample: ${file_name}\n $(date -d@${SECONDS} -u +%H:%M:%S) seconds elapsed.\n"
echo "Create haplotype file:"
echo "- bcftools query -f '[%GT ]' ${output}/${file_name}.imputed.vcf.gz | sed -E 's/[|\/]/ /g' | sed -E 's/\./?/g' | sed 's/ *$//' > ${output}/${file_name}.imputed.hap"
bcftools query -f '[%GT ]' ${output}/${file_name}.imputed.vcf.gz | sed -E 's/[|\/]/ /g' | sed -E 's/\./?/g' | sed 's/ *$//' > ${output}/${file_name}.imputed.hap
) |& tee -a "${log}"