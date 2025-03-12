#!python3

# >><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>>
# This tool sets heterozygot positions of haploid individuals to missing  >>
# Written by Demetris Taliadoros. Last update 22/03/2024                  >>
# >><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>>

# Import necessary modules

import pandas as pd
import argparse
from tqdm import tqdm

# Argument Parser
parser = argparse.ArgumentParser(description='Process heterozygous positions of haploid individuals in a VCF file.')
parser.add_argument('-v', '--vcf', type=str, help='Input VCF file', required=True)
parser.add_argument('-l', '--list', type=str, help='List of haploid individuals', required=True)
parser.add_argument('-r', '--rownum', type=int, help='Number of header rows in the VCF file minus 1', required=True)
parser.add_argument('--matt', action='store_true', help='Only output heterozygous positions without modifying the VCF')

args = parser.parse_args()

# Read the VCF file
with open(args.vcf, 'r') as vcf_file:
    header_lines = [next(vcf_file) for _ in range(args.rownum + 1)]  # Retain header
vcf = pd.read_table(args.vcf, skiprows=args.rownum)

# Read haploid individual list
with open(args.list, 'r') as file:
    haploid_individuals = set(file.read().splitlines())  # Use a set for faster lookups

# Ensure individuals exist in the VCF
haploid_individuals = [ind for ind in haploid_individuals if ind in vcf.columns]
if not haploid_individuals:
    print("Warning: No individuals in the list match the VCF column names.")

# Function to find heterozygous positions
def find_heterozygous_positions(df, individuals):
    het_positions = []
    for ind in tqdm(individuals, desc="Scanning for heterozygous positions", unit="individual"):
        for index, genotype in df[ind].items():
            if isinstance(genotype, str) and ":" in genotype:
                genotype_field = genotype.split(":")[0]  # Extract first field (genotype)
                alleles = set(genotype_field.replace("|", "/").split("/"))
                if len(alleles) > 1:  # Heterozygous check
                    het_positions.append([df.at[index, "#CHROM"], df.at[index, "POS"]])
    return pd.DataFrame(het_positions, columns=['#CHROM', 'POS']).drop_duplicates()

# If --matt is specified, output heterozygous positions without modifying the VCF
if args.matt:
    heterozygous_positions = find_heterozygous_positions(vcf, haploid_individuals)
    het_positions_file = args.vcf.replace('.vcf', '_het_positions.txt')
    heterozygous_positions.to_csv(het_positions_file, sep='\t', index=False)
    print(f"Heterozygous positions saved as {het_positions_file}")
else:
    # Function to modify heterozygous genotypes only for the affected individuals
    def set_heterozygous_to_missing(df, individuals):
        for ind in tqdm(individuals, desc="Processing individuals", unit="individual"):
            for index, genotype in df[ind].items():
                if isinstance(genotype, str) and ":" in genotype:
                    fields = genotype.split(":")  # Split genotype field
                    genotype_field = fields[0]  # Extract only genotype
                    alleles = set(genotype_field.replace("|", "/").split("/"))
                    if len(alleles) > 1:  # If heterozygous
                        fields[0] = "./."  # Replace only the genotype field
                        df.at[index, ind] = ":".join(fields)  # Reconstruct full string

    # Process the VCF
    set_heterozygous_to_missing(vcf, haploid_individuals)

    # Write the modified VCF file
    output_vcf = args.vcf.replace('.vcf', '_modified.vcf')
    with open(output_vcf, 'w') as out_vcf:
        out_vcf.writelines(header_lines)  # Preserve original header
    vcf.to_csv(output_vcf, sep='\t', index=False, mode='a')

    print(f"Modified VCF file saved as {output_vcf}")

