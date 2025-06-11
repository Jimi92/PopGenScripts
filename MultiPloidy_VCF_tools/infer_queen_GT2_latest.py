import pandas as pd
import argparse
from tqdm import tqdm
from collections import defaultdict
import gzip

# Argument Parser
parser = argparse.ArgumentParser(description='Infers queen genotype using the genotype of drone sons.')
parser.add_argument('-v', '--vcf', type=str, help='Input VCF file', required=True)
parser.add_argument('-l', '--list', type=str, help='List of individual ID and queen ID', required=True)
parser.add_argument('-r', '--rownum', type=int, help='Number of header rows in the VCF file minus 1', required=True)
parser.add_argument('-ht', '--heterozygot_thres', type=float, default=0.125, help='Threshold for heterozygosity (default=0.125)')
parser.add_argument('--haploid', action='store_true', help='Use this flag if the input VCF contains haploid genotypes')

args = parser.parse_args()




# ---------------------- PARSE INPUT FILES ---------------------------------

#  gzipped file parser
def open_maybe_gz(filename, mode='rt'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode=mode, encoding='utf-8')
    else:
        return open(filename, mode=mode, encoding='utf-8')


# Read header lines including the Sample line
with open_maybe_gz(args.vcf) as vcf_file:
    lines = vcf_file.readlines()

header_lines = lines[:args.rownum + 1]
column_line = header_lines[-1].strip().lstrip('#')
column_names = column_line.split('\t')

if len(column_names) != len(set(column_names)):
    duplicates = [x for x in column_names if column_names.count(x) > 1]
    print(f"Warning: Duplicate column names found in VCF header: {set(duplicates)}")

# Read VCF data and infers compression from filename extension
vcf = pd.read_csv(args.vcf, sep='\t', skiprows=args.rownum + 1, dtype=str, compression='infer')

# make haploid individual and queens list
droneID = []
queenID = []
with open(args.list, 'r') as file:
    for line in file:
        parts = line.strip().split()
        if len(parts) >= 2:
            droneID.append(parts[0].strip())
            queenID.append(parts[1].strip())
        else:
            print('Oops: Something is wrong with the list file. LIST should have 2 columns but have 1 for some samples')

# Check samples in VCF match the list
vcf_samples = list(vcf.columns[9:])

print("From list:", droneID)
print("From VCF :", vcf_samples)

if droneID != vcf_samples:
    print("Oops: Individuals in the list do not match the VCF column names.")
else:
    print(f"Individuals matched between list and VCF: {', '.join(droneID)}")

# ----------------- FUNCTION TO INFER QUEEN GENOTYPES ----------------------
#Open issues: 
# fix the missingness hard filter

def infer_queen_genotypes(vcf, droneID, queenID, header_lines, threshold):
    family_map = defaultdict(list)
    for drone, queen in zip(droneID, queenID):
        family_map[queen].append(drone)

    queen_genotypes = pd.DataFrame()
    queen_genotypes[vcf.columns[:9]] = vcf.iloc[:, :9]

    for queen, drones in tqdm(family_map.items(), desc="Inferring queen genotypes"):
        # Check again drones in VCF columns
        for d in drones:
            if d not in vcf.columns:
                raise ValueError(f"Drone '{d}' not found in VCF columns.")

        gt_list = []

        for _, row in vcf.iterrows():
            drone_gts = [str(row[d]).split(':')[0] for d in drones]
            # Missingness subfilter
            valid_gts = [gt for gt in drone_gts if gt not in ['./.', '.', '']]

            total = len(valid_gts)
            alt_count = sum(1 for gt in valid_gts if any(allele in gt for allele in ['1', '2']))

            if total <= 3:
                queen_gt = './.'  
            else:
                ratio = alt_count / total
                if ratio < threshold:
                    queen_gt = '0/0'
                elif ratio > (1 - threshold):
                    queen_gt = '1/1'
                else:
                    queen_gt = '0/1'

            gt_list.append(queen_gt)

        queen_genotypes[queen] = gt_list

    return queen_genotypes, header_lines

# ----------------- OUTPUT --------------------------------------------------

queen_vcf, headers = infer_queen_genotypes(vcf, droneID, queenID, header_lines, args.heterozygot_thres)

# Reconstruct the #CHROM header line with queen sample names
new_sample_header = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join(queen_vcf.columns[9:]) + '\n'
headers[-1] = new_sample_header  # Replace the old sample header line

output_vcf = args.vcf.replace('.vcf', '_queens.vcf').replace('.vcf.gz', 'GT_inferred.vcf')
with open(output_vcf, 'w') as out_vcf:
    out_vcf.writelines(headers)
    queen_vcf.to_csv(out_vcf, sep='\t', index=False, header=False)
