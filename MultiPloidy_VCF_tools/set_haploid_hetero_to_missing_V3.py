#!python3

# >><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>>
# This tool sets heterozygot positions of haploid individuals to missing  >>
# Written by Demetris Taliadoros. Last update 19/01/2026                  >>
# >><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>>

#!/usr/bin/env python3
"""
This tool sets heterozygous positions of haploid individuals to missing.
Written by Demetris Taliadoros. Last update 22/03/2024
Updated: multi-allelic AD support + bgzip/gzip VCF reading + robust AD parsing via FORMAT.
"""

import pandas as pd
import argparse
import re
import gzip
from tqdm import tqdm
from typing import Optional, List


# -------------------------
# Helpers for (b)gzip VCF I/O
# -------------------------
def open_maybe_gzip(path: str, mode: str = "rt"):
    """
    Open plain text or gzip/bgzip-compressed files transparently.
    bgzip files are gzip-compatible at the stream level, so gzip.open works.
    """
    lower = path.lower()
    if lower.endswith((".gz", ".bgz", ".bgzip")):
        return gzip.open(path, mode)
    return open(path, mode)


def pandas_read_vcf_table(path: str, skiprows: int) -> pd.DataFrame:
    """
    Read VCF body with pandas, handling gzip/bgzip compression.
    """
    lower = path.lower()
    if lower.endswith((".gz", ".bgz", ".bgzip")):
        return pd.read_table(path, skiprows=skiprows, compression="gzip")
    return pd.read_table(path, skiprows=skiprows)


# -------------------------
# Genotype checks
# -------------------------
def is_het_from_gt(cell: str) -> bool:
    """
    GT-based heterozygosity detection:
    - Looks at first FORMAT subfield (GT).
    - Splits by / or |.
    - If >1 distinct allele codes => heterozygous (e.g., 0/1, 0/2, 1/2).
    """
    if not isinstance(cell, str) or ":" not in cell:
        return False
    gt = cell.split(":", 1)[0]
    alleles = set(re.split(r"[|/]", gt))
    # ignore missing tokens
    alleles.discard(".")
    return len(alleles) > 1


def parse_ad_from_row(format_str: str, sample_str: str) -> Optional[List[int]]:
    """
    Parse AD from a VCF row using FORMAT to locate AD field.
    Returns list of ints: [ref_depth, alt1_depth, alt2_depth, ...]
    """
    if not isinstance(format_str, str) or not isinstance(sample_str, str):
        return None
    if ":" not in sample_str:
        return None

    fmt_keys = format_str.split(":")
    vals = sample_str.split(":")
    if len(vals) < len(fmt_keys):
        # Some VCFs truncate trailing missing fields; handle safely by padding.
        vals = vals + [""] * (len(fmt_keys) - len(vals))

    try:
        ad_idx = fmt_keys.index("AD")
    except ValueError:
        return None  # AD not present in FORMAT

    ad_raw = vals[ad_idx]
    if not ad_raw or ad_raw == ".":
        return None

    parts = ad_raw.split(",")
    depths: List[int] = []
    for p in parts:
        if p in ("", "."):
            return None
        try:
            depths.append(int(p))
        except ValueError:
            return None

    if len(depths) < 2:
        # AD should usually have at least ref + one alt
        return None
    return depths


def ad_has_alt_balanced_against_ref(ad_depths: List[int], low: float = 0.2, high: float = 1.8) -> bool:
    """
    Multi-allelic aware rule:
    Compare each ALT depth against REF depth.

    - ad_depths = [ref, alt1, alt2, ...]
    - If ref > 0: compute alt/ref ratio for each alt; if any ratio in [low, high] => flag
    - If ref == 0: if any alt > 0 => flag (since ratio undefined but evidence of non-ref)
    """
    ref = ad_depths[0]
    alts = ad_depths[1:]

    if ref <= 0:
        return any(a > 0 for a in alts)

    for alt in alts:
        if alt <= 0:
            continue
        ratio = alt / ref
        if low <= ratio <= high:
            return True
    return False


# -------------------------
# Main scanning / modifying
# -------------------------
def find_heterozygous_positions_gt(df: pd.DataFrame, individuals: List[str]) -> pd.DataFrame:
    het_positions = []
    for ind in tqdm(individuals, desc="Scanning for heterozygous positions (GT)", unit="individual"):
        for idx, cell in df[ind].items():
            if is_het_from_gt(cell):
                het_positions.append([df.at[idx, "#CHROM"], df.at[idx, "POS"]])
    return pd.DataFrame(het_positions, columns=["#CHROM", "POS"]).drop_duplicates()


def find_positions_ad(df: pd.DataFrame, individuals: List[str], low: float = 0.2, high: float = 1.8) -> pd.DataFrame:
    ad_positions = []
    # We use the row's FORMAT column to locate AD for each sample.
    if "FORMAT" not in df.columns:
        raise ValueError("VCF body is missing FORMAT column; cannot use --AD mode.")

    for ind in tqdm(individuals, desc="Scanning positions (AD vs REF)", unit="individual"):
        for idx, sample_cell in df[ind].items():
            fmt = df.at[idx, "FORMAT"]
            ad = parse_ad_from_row(fmt, sample_cell)
            if ad is None:
                continue
            if ad_has_alt_balanced_against_ref(ad, low=low, high=high):
                ad_positions.append([df.at[idx, "#CHROM"], df.at[idx, "POS"]])

    return pd.DataFrame(ad_positions, columns=["#CHROM", "POS"]).drop_duplicates()


def set_positions_to_missing(df: pd.DataFrame, individuals: List[str], use_ad: bool, low: float = 0.2, high: float = 1.8) -> None:
    """
    Set flagged haploid sample genotypes to missing:
    - If use_ad: flag rows where AD indicates alt depth balanced vs ref depth (multi-allelic aware)
    - Else: flag rows where GT is heterozygous (contains >1 distinct allele)
    """
    if use_ad and "FORMAT" not in df.columns:
        raise ValueError("VCF body is missing FORMAT column; cannot use --AD mode.")

    for ind in tqdm(individuals, desc="Processing individuals", unit="individual"):
        for idx, cell in df[ind].items():
            if not isinstance(cell, str) or ":" not in cell:
                continue

            flag = False
            if use_ad:
                fmt = df.at[idx, "FORMAT"]
                ad = parse_ad_from_row(fmt, cell)
                if ad is not None and ad_has_alt_balanced_against_ref(ad, low=low, high=high):
                    flag = True
            else:
                flag = is_het_from_gt(cell)

            if flag:
                fields = cell.split(":")
                fields[0] = "./."
                df.at[idx, ind] = ":".join(fields)


def main():
    parser = argparse.ArgumentParser(
        description="Process heterozygous positions of haploid individuals in a VCF file."
    )
    parser.add_argument("-v", "--vcf", type=str, required=True, help="Input VCF file (.vcf, .vcf.gz, .vcf.bgz)")
    parser.add_argument("-l", "--list", type=str, required=True, help="List of haploid individuals")
    parser.add_argument("-r", "--rownum", type=int, required=True, help="Number of header rows in the VCF file minus 1")
    parser.add_argument("--matt", action="store_true", help="Only output positions without modifying the VCF")
    parser.add_argument("--AD", action="store_true", help="Use AD field to check allele balance instead of GT field")
    parser.add_argument("--low", type=float, default=0.2, help="Lower bound for ALT/REF ratio in AD mode (default 0.2)")
    parser.add_argument("--high", type=float, default=1.8, help="Upper bound for ALT/REF ratio in AD mode (default 1.8)")

    args = parser.parse_args()

    # Read header lines (works for plain and gz/bgzip)
    with open_maybe_gzip(args.vcf, "rt") as vcf_file:
        header_lines = [next(vcf_file) for _ in range(args.rownum + 1)]

    # Read body table
    vcf = pandas_read_vcf_table(args.vcf, skiprows=args.rownum)

    # Read haploid individual list
    with open(args.list, "r") as f:
        haploid_individuals = set(f.read().splitlines())

    # Ensure individuals exist in the VCF
    haploid_individuals = [ind for ind in haploid_individuals if ind in vcf.columns]
    if not haploid_individuals:
        print("Error: No valid haploid individuals found in VCF. Exiting.")
        raise SystemExit(1)

    # If --matt specified: output positions only
    if args.matt:
        if args.AD:
            positions = find_positions_ad(vcf, haploid_individuals, low=args.low, high=args.high)
            out = re.sub(r"(\.vcf)(\.(gz|bgz|bgzip))?$", r"_AD_positions.txt", args.vcf, flags=re.IGNORECASE)
            positions.to_csv(out, sep="\t", index=False)
            print(f"AD-based positions saved as {out}")
        else:
            positions = find_heterozygous_positions_gt(vcf, haploid_individuals)
            out = re.sub(r"(\.vcf)(\.(gz|bgz|bgzip))?$", r"_het_positions.txt", args.vcf, flags=re.IGNORECASE)
            positions.to_csv(out, sep="\t", index=False)
            print(f"Heterozygous positions saved as {out}")
        raise SystemExit(0)

    # Modify and write VCF (output is uncompressed .vcf)
    set_positions_to_missing(vcf, haploid_individuals, use_ad=args.AD, low=args.low, high=args.high)

    output_vcf = re.sub(r"(\.vcf)(\.(gz|bgz|bgzip))?$", r"_modified.vcf", args.vcf, flags=re.IGNORECASE)
    with open(output_vcf, "w") as out_vcf:
        out_vcf.writelines(header_lines)
    vcf.to_csv(output_vcf, sep="\t", index=False, mode="a", header=True)

    print(f"Modified VCF file saved as {output_vcf}")


if __name__ == "__main__":
    main()

