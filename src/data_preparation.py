import pandas as pd
from utils import get_data
import argparse

parser = argparse.ArgumentParser(
                    prog='GIP',
                    description='Data generation',
                    epilog='Help')

parser.add_argument('data_name', metavar='n', type=str, help='Prefix of data name')
parser.add_argument('rate', metavar='r', type=str, help='Missing rate')
args = parser.parse_args()

rate = float(args.rate)
vcf_path = "../data/HLA/" + args.data_name + ".vcf"
output_VCF = "../data/HLA/" + args.data_name + "." + str(rate) + ".missing.vcf"
output_mask = "../data/HLA/" + args.data_name + "." + str(rate) + ".mask.csv"


# Read a vcf file
def read_vcf(vcf_path):
    with open(vcf_path, "rt") as ifile:
          for line in ifile:
            if line.startswith("#CHROM"):
                  vcf_names = [x for x in line.split('\t')]
                  break
    ifile.close()
    data = pd.read_csv(vcf_path, comment='#', delim_whitespace=True, header=None, names=vcf_names)
    return data

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
    geno = data.iloc[:, 9::]
    geno, geno_miss, m = get_data(geno, rate)
    gn = pd.concat([data.iloc[:, 0:9], geno_miss], axis = 1)
    gn.rename(columns={'NA21144\n':'NA21144'}, inplace=True)
    
    # Header for vcf file
    header = """##fileformat=VCFv4.1
##FILTER=<ID=PASS,Description="All filters passed">
##fileDate=20150218
##reference=ftp://ftp.1000genomes.ebi.ac.uk//vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
##source=1000GenomesPhase3Pipeline
##contig=<ID=chr1,assembly=b37,length=249250621>
##contig=<ID=chr2,assembly=b37,length=243199373>
##contig=<ID=chr3,assembly=b37,length=198022430>
##contig=<ID=chr4,assembly=b37,length=191154276>
##contig=<ID=chr5,assembly=b37,length=180915260>
##contig=<ID=chr6,assembly=b37,length=171115067>
##contig=<ID=chr7,assembly=b37,length=159138663>
##contig=<ID=chr8,assembly=b37,length=146364022>
##contig=<ID=chr9,assembly=b37,length=141213431>
##contig=<ID=chr10,assembly=b37,length=135534747>
##contig=<ID=chr11,assembly=b37,length=135006516>
##contig=<ID=chr12,assembly=b37,length=133851895>
##contig=<ID=chr13,assembly=b37,length=115169878>
##contig=<ID=chr14,assembly=b37,length=107349540>
##contig=<ID=chr15,assembly=b37,length=102531392>
##contig=<ID=chr16,assembly=b37,length=90354753>
##contig=<ID=chr17,assembly=b37,length=81195210>
##contig=<ID=chr18,assembly=b37,length=78077248>
##contig=<ID=chr19,assembly=b37,length=59128983>
##contig=<ID=chr20,assembly=b37,length=63025520>
##contig=<ID=chr21,assembly=b37,length=48129895>
##contig=<ID=chr22,assembly=b37,length=51304566>
##contig=<ID=GL000191.1,assembly=b37,length=106433>
##contig=<ID=GL000192.1,assembly=b37,length=547496>
##contig=<ID=GL000193.1,assembly=b37,length=189789>
##contig=<ID=GL000194.1,assembly=b37,length=191469>
##contig=<ID=GL000195.1,assembly=b37,length=182896>
##contig=<ID=GL000196.1,assembly=b37,length=38914>
##contig=<ID=GL000197.1,assembly=b37,length=37175>
##contig=<ID=GL000198.1,assembly=b37,length=90085>
##contig=<ID=GL000199.1,assembly=b37,length=169874>
##contig=<ID=GL000200.1,assembly=b37,length=187035>
##contig=<ID=GL000201.1,assembly=b37,length=36148>
##contig=<ID=GL000202.1,assembly=b37,length=40103>
##contig=<ID=GL000203.1,assembly=b37,length=37498>
##contig=<ID=GL000204.1,assembly=b37,length=81310>
##contig=<ID=GL000205.1,assembly=b37,length=174588>
##contig=<ID=GL000206.1,assembly=b37,length=41001>
##contig=<ID=GL000207.1,assembly=b37,length=4262>
##contig=<ID=GL000208.1,assembly=b37,length=92689>
##contig=<ID=GL000209.1,assembly=b37,length=159169>
##contig=<ID=GL000210.1,assembly=b37,length=27682>
##contig=<ID=GL000211.1,assembly=b37,length=166566>
##contig=<ID=GL000212.1,assembly=b37,length=186858>
##contig=<ID=GL000213.1,assembly=b37,length=164239>
##contig=<ID=GL000214.1,assembly=b37,length=137718>
##contig=<ID=GL000215.1,assembly=b37,length=172545>
##contig=<ID=GL000216.1,assembly=b37,length=172294>
##contig=<ID=GL000217.1,assembly=b37,length=172149>
##contig=<ID=GL000218.1,assembly=b37,length=161147>
##contig=<ID=GL000219.1,assembly=b37,length=179198>
##contig=<ID=GL000220.1,assembly=b37,length=161802>
##contig=<ID=GL000221.1,assembly=b37,length=155397>
##contig=<ID=GL000222.1,assembly=b37,length=186861>
##contig=<ID=GL000223.1,assembly=b37,length=180455>
##contig=<ID=GL000224.1,assembly=b37,length=179693>
##contig=<ID=GL000225.1,assembly=b37,length=211173>
##contig=<ID=GL000226.1,assembly=b37,length=15008>
##contig=<ID=GL000227.1,assembly=b37,length=128374>
##contig=<ID=GL000228.1,assembly=b37,length=129120>
##contig=<ID=GL000229.1,assembly=b37,length=19913>
##contig=<ID=GL000230.1,assembly=b37,length=43691>
##contig=<ID=GL000231.1,assembly=b37,length=27386>
##contig=<ID=GL000232.1,assembly=b37,length=40652>
##contig=<ID=GL000233.1,assembly=b37,length=45941>
##contig=<ID=GL000234.1,assembly=b37,length=40531>
##contig=<ID=GL000235.1,assembly=b37,length=34474>
##contig=<ID=GL000236.1,assembly=b37,length=41934>
##contig=<ID=GL000237.1,assembly=b37,length=45867>
##contig=<ID=GL000238.1,assembly=b37,length=39939>
##contig=<ID=GL000239.1,assembly=b37,length=33824>
##contig=<ID=GL000240.1,assembly=b37,length=41933>
##contig=<ID=GL000241.1,assembly=b37,length=42152>
##contig=<ID=GL000242.1,assembly=b37,length=43523>
##contig=<ID=GL000243.1,assembly=b37,length=43341>
##contig=<ID=GL000244.1,assembly=b37,length=39929>
##contig=<ID=GL000245.1,assembly=b37,length=36651>
##contig=<ID=GL000246.1,assembly=b37,length=38154>
##contig=<ID=GL000247.1,assembly=b37,length=36422>
##contig=<ID=GL000248.1,assembly=b37,length=39786>
##contig=<ID=GL000249.1,assembly=b37,length=38502>
##contig=<ID=MT,assembly=b37,length=16569>
##contig=<ID=NC_007605,assembly=b37,length=171823>
##contig=<ID=chrX,assembly=b37,length=155270560>
##contig=<ID=Y,assembly=b37,length=59373566>
##contig=<ID=hs37d5,assembly=b37,length=35477943>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">
##INFO=<ID=CS,Number=1,Type=String,Description="Source call set.">
##INFO=<ID=END,Number=1,Type=Integer,Description="End coordinate of this variant">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
##INFO=<ID=MC,Number=.,Type=String,Description="Merged calls.">
##INFO=<ID=MEINFO,Number=4,Type=String,Description="Mobile element info of the form NAME,START,END<POLARITY; If there is only 5' OR 3' support for this call, will be NULL NULL for START and END">
##INFO=<ID=MEND,Number=1,Type=Integer,Description="Mitochondrial end coordinate of inserted sequence">
##INFO=<ID=MLEN,Number=1,Type=Integer,Description="Estimated length of mitochondrial insert">
##INFO=<ID=MSTART,Number=1,Type=Integer,Description="Mitochondrial start coordinate of inserted sequence">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="SV length. It is only calculated for structural variation MEIs. For other types of SVs; one may calculate the SV length by INFO:END-START+1, or by finding the difference between lengthes of REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=TSD,Number=1,Type=String,Description="Precise Target Site Duplication for bases, if unknown, value will be NULL">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Total number of alternate alleles in called genotypes">
##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated allele frequency in the range (0,1)">
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=EAS_AF,Number=A,Type=Float,Description="Allele frequency in the EAS populations calculated from AC and AN, in the range (0,1)">
##INFO=<ID=EUR_AF,Number=A,Type=Float,Description="Allele frequency in the EUR populations calculated from AC and AN, in the range (0,1)">
##INFO=<ID=AFR_AF,Number=A,Type=Float,Description="Allele frequency in the AFR populations calculated from AC and AN, in the range (0,1)">
##INFO=<ID=AMR_AF,Number=A,Type=Float,Description="Allele frequency in the AMR populations calculated from AC and AN, in the range (0,1)">
##INFO=<ID=SAS_AF,Number=A,Type=Float,Description="Allele frequency in the SAS populations calculated from AC and AN, in the range (0,1)">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth; only low coverage data were counted towards the DP, exome data were not used">
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele. Format: AA|REF|ALT|IndelType. AA: Ancestral allele, REF:Reference Allele, ALT:Alternate Allele, IndelType:Type of Indel (REF, ALT and IndelType are only defined for indels)">
##INFO=<ID=VT,Number=.,Type=String,Description="indicates what type of variant the line represents">
##INFO=<ID=EX_TARGET,Number=0,Type=Flag,Description="indicates whether a variant is within the exon pull down target boundaries">
##INFO=<ID=MULTI_ALLELIC,Number=0,Type=Flag,Description="indicates whether a site is multi-allelic">
##INFO=<ID=STRAND_FLIP,Number=0,Type=Flag,Description="Indicates that the reference strand has changed between GRCh37 and GRCh38">
##INFO=<ID=REF_SWITCH,Number=0,Type=Flag,Description="Indicates that the reference allele has changed">
##INFO=<ID=DEPRECATED_RSID,Number=.,Type=String,Description="dbsnp rs IDs that have been merged into other rs IDs or do not map to GRCh38">
##INFO=<ID=RSID_REMOVED,Number=.,Type=String,Description="dbsnp rs IDs removed from this variant, due to either the variant splitting up or being deprecated/merged">
##INFO=<ID=GRCH37_38_REF_STRING_MATCH,Number=0,Type=Flag,Description="Indicates reference allele in origin GRCh37 vcf string-matches reference allele in dbsnp GRCh38 vcf">
##INFO=<ID=NOT_ALL_RSIDS_STRAND_CHANGE_OR_REF_SWITCH,Number=0,Type=Flag,Description="Indicates only some of the rs IDs in origin GRCh37 vcf switched strands or switched strands and changed reference allele. This would result in rs IDs being split into multiple lines">
##INFO=<ID=GRCH37_POS,Number=1,Type=Integer,Description="Position in origin GRCh37 vcf">
##INFO=<ID=GRCH37_REF,Number=1,Type=String,Description="Representation of reference allele in origin GRCh37 vcf">
##INFO=<ID=ALLELE_TRANSFORM,Number=0,Type=Flag,Description="Indicates that at least some of the alleles have changed in how they're represented, e.g. through left shifting.">
##INFO=<ID=REF_NEW_ALLELE,Number=0,Type=Flag,Description="Indicates that the reference allele is an allele not present in the origin GRCh37 vcf">
##INFO=<ID=CHROM_CHANGE_BETWEEN_ASSEMBLIES,Number=.,Type=String,Description="dbsnp rs IDs that are mapped to a different chromosome between GRCh37 and GRCh38">
##bcftools_annotateVersion=1.9+htslib-1.9
##bcftools_annotateCommand=annotate --rename-chrs chr_names.txt -Ou ALL.chr6_GRCh38.genotypes.20170504.vcf.gz; Date=Mon Aug 28 20:38:30 2023
##bcftools_viewVersion=1.9+htslib-1.9
##bcftools_viewCommand=view -e 'INFO/AC<3 | INFO/AN-INFO/AC<3' -Ou; Date=Mon Aug 28 20:38:30 2023
##bcftools_normVersion=1.9+htslib-1.9
##bcftools_normCommand=norm -m -any -Ou; Date=Mon Aug 28 20:38:30 2023
##bcftools_viewCommand=view -i 'INFO/VT="SNP" | INFO/VT="INDEL"' -Ou; Date=Mon Aug 28 20:38:30 2023
##bcftools_normCommand=norm -f ../ref_genome/hg38.fasta -d none -Ou; Date=Mon Aug 28 20:38:30 2023
##bcftools_viewCommand=view -m 2 -M 2 -Ou; Date=Mon Aug 28 20:38:30 2023
##bcftools_viewCommand=view -g ^miss -Oz -o 1000GP_chr6.vcf.gz; Date=Mon Aug 28 20:38:30 2023
##bcftools_viewCommand=view -e ID=@1000GP_chr6.dup_id -Oz -o 1000GP_filtered_chr6.vcf.gz 1000GP_chr6.vcf.gz; Date=Tue Aug 29 18:58:01 2023
##bcftools_pluginVersion=1.9+htslib-1.9
##bcftools_pluginCommand=plugin fill-tags -Oz -o 1000GP_AF_chr6.vcf.gz -- 1000GP_filtered_chr6.vcf.gz -t AF; Date=Sun Sep 24 18:38:53 2023
##bcftools_viewVersion=1.14+htslib-1.15.1
##bcftools_viewCommand=view -R HLA1_chr6.bed 1000GP_AF_chr6.vcf.gz; Date=Thu Mar 28 08:24:00 2024
"""
    print("********************************** SAVING **********************************")
    with open(output_VCF, 'w') as vcf:
        vcf.write(header)
    gn.to_csv(output_VCF, sep="\t", mode='a', index=False)
    pd.DataFrame(columns = list(gn.columns)[9::], data = m).to_csv(output_mask)
    print("")
    print("********************************* FINISHED *********************************")
    print("")
    