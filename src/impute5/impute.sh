#!/bin/bash
ref_panel=$1 # /work/users/minhnth/obesity/imputation/ref_panel_for_minimac/1000GP_AF_chr6.msav
input_dir=$2 # ../data/HLA_3_alleles/output/HLA1_chr6.SAS.0.25.missing.imputed.minimac.vcf
file_name=$3 # prefix of filename
output_dir=$4 # ../data/obesity/imputation_output/
region=$5 # chr6:292100-351355
genetic_map=$6
num_threads=$7
impute5_path=$8 #GIP/impute5_v1.2.0/impute5_v1.2.0_static

echo "Waiting for Impute5"
rm ${output_dir}/${file_name}.imputed.ip5.hap
: << 'COMMENT'
# fetch shapeit4_map data for impute5
shapeit4_map_dir=imputation/genetic_map_for_imputation/
mkdir -p ${shapeit4_map_dir}
cd ${shapeit4_map_dir}
wget -q -O- https://github.com/odelaneau/shapeit4/raw/master/maps/genetic_maps.b38.tar.gz | tar xzvf -
COMMENT

mkdir -p ${output_dir}
#cp ../../data/HLA/*.vcf impute5_test
echo "" > ${file_name}.log.txt
log=$(readlink -f ${file_name}.log.txt)
(
#cd impute5_test

#for i in *.missing.vcf ; do
echo "Sample: ${file_name}"
cat ${input_dir}/${file_name}.vcf | sed -r 's:\.\|\.:0\|0:g' > ${input_dir}/${file_name}.ip5.vcf # Change .|. to 0|0
bcftools view ${input_dir}/${file_name}.ip5.vcf -O z -o ${input_dir}/${file_name}.ip5.vcf.gz
tabix -f ${input_dir}/${file_name}.ip5.vcf.gz

#cat $i | sed -r 's:chr::g' > $i.new.vcf # remove chr
#mv $i.new.vcf $i
echo "${impute5_path} --h ${ref_panel} --m ${genetic_map} --g ${input_dir}/${file_name}.ip5.vcf.gz --r ${region} --o ${input_dir}/${file_name}.imputed.ip5.vcf --buffer-region ${region}" >> $log
# ${impute5_path} --h ${reference_path} --m ${shapeit4_map} --g $i.gz --r 6:292100-351355 --o $i.imputed.vcf --buffer-region 6:0-400000 >> $log

SECONDS=0

${impute5_path} --h ${ref_panel} \
--m ${genetic_map} \
--g ${input_dir}/${file_name}.ip5.vcf.gz \
--r ${region} \
--buffer-region ${region} \
--o ${output_dir}/${file_name}.imputed.ip5.info.vcf

duration=$SECONDS
echo "$((duration / 60)) minutes and $((duration % 60)) seconds elapsed."

echo "Sample: ${file_name}" >> $log
echo "$((duration / 60)) minutes and $((duration % 60)) seconds elapsed." >> $log
bcftools annotate -x ^FORMAT/GT ${output_dir}/${file_name}.imputed.ip5.info.vcf > ${output_dir}/${file_name}.imputed.ip5.vcf
#done
)
#echo "Evaluate: "
#python impute_evaluate.py >> $log
echo "Saving ..."
# Get haploid data from Impute5's result
python gip/vcf2hap.py \
${input_dir}/${file_name:0:-7}ori.vcf \
${input_dir}/${file_name}.vcf \
${output_dir}/${file_name}.imputed.ip5.vcf \
${output_dir}/${file_name}.imputed.ip5.hap \
impute5

echo "Completed"