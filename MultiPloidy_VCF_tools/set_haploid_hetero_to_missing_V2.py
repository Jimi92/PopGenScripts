#!python3

# >><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>>
# This tool sets heterozygot positions of haploid individuals to missing  >>
# Written by Demetris Taliadoros. Last update 22/03/2024                  >>
# >><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>>

#!python3

import pandas as pd
import argparse
import re
from tqdm import tqdm

# Argument Parser
parser = argparse.ArgumentParser(description='Process heterozygous positions of haploid individuals in a VCF file.')
parser.add_argument('-v', '--vcf', type=str, help='Input VCF file', required=True)
parser.add_argument('-l', '--list', type=str, help='List of haploid individuals', required=True)
parser.add_argument('-r', '--rownum', type=int, help='Number of header rows in the VCF file minus 1', required=True)
parser.add_argument('--matt', action='store_true', help='Only output positions without modifying the VCF')
parser.add_argument('--AD', action='store_true', help='Use AD field to check allele balance instead of GT field')

args = parser.parse_args()

# Read the VCF file
with open(args.vcf, 'r') as vcf_file:
    header_lines = [next(vcf_file) for _ in range(args.rownum + 1)]  # Retain header
vcf = pd.read_table(args.vcf, skiprows=args.rownum)

# Read haploid individual list
with open(args.list, 'r') as file:
    haploid_individuals = set(file.read().splitlines())

# Ensure individuals exist in the VCF
haploid_individuals = [ind for ind in haploid_individuals if ind in vcf.columns]
if not haploid_individuals:
    print("Error: No valid haploid individuals found in VCF. Exiting.")
    exit(1)

# Function to check heterozygosity via GT field
def find_heterozygous_positions(df, individuals):
    het_positions = []
    for ind in tqdm(individuals, desc="Scanning for heterozygous positions", unit="individual"):
        for index, genotype in df[ind].items():
            if isinstance(genotype, str) and ":" in genotype:
                genotype_field = genotype.split(":")[0]
                alleles = set(re.split(r'[|/]', genotype_field))
                if len(alleles) > 1:  # Heterozygous check
                    het_positions.append([df.at[index, "#CHROM"], df.at[index, "POS"]])
    return pd.DataFrame(het_positions, columns=['#CHROM', 'POS']).drop_duplicates()

# Function to check allele balance via AD field
def find_AD_outliers(df, individuals):
    ad_positions = []
    for ind in tqdm(individuals, desc="Checking AD values", unit="individual"):
        for index, genotype in df[ind].items():
            if isinstance(genotype, str) and ":" in genotype:
                fields = genotype.split(":")
                if len(fields) > 1:  # Ensure AD field exists
                    ad_values = fields[1]  # AD field (second field)
                    if "," in ad_values:
                        try:
                            AD1, AD2 = map(int, ad_values.split(","))
                            if AD1 > 0:  # Avoid division by zero
                                ratio = AD2 / AD1
                                if 0.2 <= ratio <= 1.8:  # Check ratio
                                    ad_positions.append([df.at[index, "#CHROM"], df.at[index, "POS"]])
                        except ValueError:
                            continue  # Skip invalid entries
    return pd.DataFrame(ad_positions, columns=['#CHROM', 'POS']).drop_duplicates()

# If --matt is specified, output the positions only
if args.matt:
    if args.AD:
        ad_positions = find_AD_outliers(vcf, haploid_individuals)
        ad_file = args.vcf.replace('.vcf', '_AD_positions.txt')
        ad_positions.to_csv(ad_file, sep='\t', index=False)
        print(f"AD-based positions saved as {ad_file}")
    else:
        heterozygous_positions = find_heterozygous_positions(vcf, haploid_individuals)
        het_positions_file = args.vcf.replace('.vcf', '_het_positions.txt')
        heterozygous_positions.to_csv(het_positions_file, sep='\t', index=False)
        print(f"Heterozygous positions saved as {het_positions_file}")
    exit(0)  # Exit after saving the positions

# Function to modify heterozygous genotypes
def set_heterozygous_to_missing(df, individuals):
    for ind in tqdm(individuals, desc="Processing individuals", unit="individual"):
        for index, genotype in df[ind].items():
            if isinstance(genotype, str) and ":" in genotype:
                fields = genotype.split(":")
                genotype_field = fields[0]
                alleles = set(re.split(r'[|/]', genotype_field))
                if len(alleles) > 1:
                    fields[0] = "./."
                    df.at[index, ind] = ":".join(fields)

# Modify heterozygous genotypes and save VCF
set_heterozygous_to_missing(vcf, haploid_individuals)

output_vcf = args.vcf.replace('.vcf', '_modified.vcf')
with open(output_vcf, 'w') as out_vcf:
    out_vcf.writelines(header_lines)
vcf.to_csv(output_vcf, sep='\t', index=False, mode='a', header=True)

print(f"Modified VCF file saved as {output_vcf}")
