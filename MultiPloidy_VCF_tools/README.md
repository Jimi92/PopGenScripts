# MultiPloidy_VCF_tools

These scripts were designed to deal with vcf that contain individuals with different ploidies!



## set_hetero_haploid_to_missing

This tool sets the heterozygous positions of haploid indivisuals to missing or prints a list with these positions that can be used for masking.

Necessary inputs:

 * A VCF file

 * A list of the haploid individuals in the vcf file. One individual per line.
  
 * The number of header lines in the vcf file.
  
## to set heterozygous positions to missing

Usage: set_haploid_hetero_to_missing.py -vcf file.vcf -l list_of_haploids_in_vcf.txt -r number_of_vcf_heder_lines_minus_1


## to get a list of heterozygous positions add the flag [--matt]

Usage: set_haploid_hetero_to_missing.py -vcf file.vcf -l list_of_haploids_in_vcf.txt -r number_of_vcf_heder_lines_minus_1 --matt


Comments:

This script can handle multiallelic positions.
