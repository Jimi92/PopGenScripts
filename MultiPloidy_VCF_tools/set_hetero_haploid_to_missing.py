##############################################################
##                                                          ##
## Written by Demetris Taliadoros. Last updated 22/03/2024  ##
##                                                          ##
##############################################################

# Import necessary modules

import pandas as pd
import numpy as np
import argparse
from tqdm import tqdm

## >><<>><<>><<>><<>><<>><<>><<

# Create parser, Add arguments 

parser = argparse.ArgumentParser(description='Set male heterozygot positions to missing')
parser.add_argument('-v','--vcf', type=str, help='Input VCF', required=True)
parser.add_argument('-l','--list',type=str, help='List of haploids to be processed', required=True)
parser.add_argument('-r','--rownum', type=int, help='Number of header rows in the vcf file minus 1', required=True)
args = parser.parse_args()

# >><<>><<>><<>><<>><<>><<>><<

# Main body

# Read VCF
vcf = pd.read_table(args.vcf, skiprows=args.rownum)

# Read list of haploids
with open(args.list, 'r') as file:
    list_of_males = file.read().splitlines()

# Function to set heterozygous positions to missing for every haploid individual
def filter_homozygous(df,list_of_males):
    for c in tqdm(list_of_males, desc="Processing haploids", unit="haploid"):
        df[c].replace('0/1','./.', inplace=True)
        df[c].replace('1/0','./.', inplace=True)
        df[c].replace('0|1','./.', inplace=True)
        df[c].replace('1|0','./.', inplace=True)

# Call the function
filter_homozygous(vcf, list_of_males)

# Write new file
modified_vcf_file = args.vcf.replace('.vcf', '_modified.vcf')
vcf.to_csv(modified_vcf_file, sep='\t', index=False, header=False)

print(f"Modified VCF file saved as {modified_vcf_file}")
