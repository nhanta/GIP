impute5_path=/home/pbthang/GIP/impute5_v1.2.0/impute5_v1.2.0_static
shapeit4_map=/home/pbthang/GIP/imputation/genetic_map_for_imputation/shapeit4/chr6.b38.gmap.gz
reference_path=/home/pbthang/GIP/imputation/ref_panel/ALL.chr6_GRCh38.genotypes.20170504.vcf.gz


: << 'COMMENT'
# fetch shapeit4_map data for impute5
shapeit4_map_dir=imputation/genetic_map_for_imputation/
mkdir -p ${shapeit4_map_dir}
cd ${shapeit4_map_dir}
wget -q -O- https://github.com/odelaneau/shapeit4/raw/master/maps/genetic_maps.b38.tar.gz | tar xzvf -
COMMENT

mkdir -p impute5_test
cp ../../data/HLA/*.vcf impute5_test
echo "" > log.txt
log=$(readlink -f log.txt)
(
cd impute5_test
for i in *.missing.vcf ; do
  echo "Sample: $i"
  cat $i | sed -r 's:chr::g' | sed -r 's:\.\|\.:0\|0:g' > $i.new.vcf # remove chr and dont know why i must remove .|.
  bcftools view $i.new.vcf -O z -o $i.gz
  tabix -f $i.gz
  cat $i | sed -r 's:chr::g' > $i.new.vcf # remove chr
  mv $i.new.vcf $i
  echo "${impute5_path} --h ${reference_path} --m ${shapeit4_map} --g $i.gz --r 6:292100-351355 --o $i.imputed.vcf --buffer-region 6:0-400000" >> $log
  ${impute5_path} --h ${reference_path} --m ${shapeit4_map} --g $i.gz --r 6:292100-351355 --o $i.imputed.vcf --buffer-region 6:0-400000 >> $log
  echo "" >> $log
  bcftools annotate -x ^FORMAT/GT $i.imputed.vcf > $i.final.vcf
done
)
echo "Evaluate: "
python impute_evaluate.py >> $log

