bcftools query -f '[%GT ]' $1 \
  | sed -E 's/[|\/]/ /g' | sed -E 's/\./?/g' | sed 's/ *$//' \
  > $1.hap