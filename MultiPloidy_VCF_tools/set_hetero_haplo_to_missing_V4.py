#!/usr/bin/env python3

# >><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>>
# This tool sets heterozygous positions of haploid individuals to missing.
# Optimized streaming version:
#   - No pandas
#   - Low-memory streaming
#   - GT and AD modes
#   - Multi-allelic AD support
#   - Ignores AD positions where ANY AD value is zero
#   - Reads plain VCF / gzip / bgzip
#   - Writes proper multithreaded BGZF VCF output
# >><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>>

import argparse
import gzip
import os
import re
import shutil
import subprocess
import sys
from typing import Optional, List, Dict, TextIO


# =============================================================================
# File handling
# =============================================================================

def open_maybe_gzip(path: str, mode: str = "rt"):
    """
    Transparently open plain-text, gzip, or bgzip-compressed input.

    BGZF is gzip-compatible for sequential reading, therefore gzip.open()
    is sufficient when random access is not required.
    """
    lower = path.lower()

    if lower.endswith((".gz", ".bgz", ".bgzip")):
        return gzip.open(path, mode)

    return open(path, mode)


def make_output_filename(input_vcf: str, suffix: str) -> str:
    """
    Remove .vcf / .vcf.gz / .vcf.bgz / .vcf.bgzip and append suffix.
    """
    base = re.sub(
        r"\.vcf(?:\.(?:gz|bgz|bgzip))?$",
        "",
        input_vcf,
        flags=re.IGNORECASE,
    )

    # If the filename didn't match the usual VCF endings, still create
    # a sensible output filename.
    if base == input_vcf:
        base = input_vcf

    return base + suffix


class BgzipWriter:
    """
    Context manager for writing proper BGZF output via the external
    bgzip executable.

    Text is written to .write(), encoded to bytes, then sent to bgzip.
    """

    def __init__(self, path: str, threads: int = 4):
        self.path = path
        self.threads = threads
        self.output_handle = None
        self.process = None

    def __enter__(self):
        if shutil.which("bgzip") is None:
            raise RuntimeError(
                "bgzip was not found in PATH.\n"
                "Install HTSlib/tabix/bcftools and ensure 'bgzip' is available."
            )

        self.output_handle = open(self.path, "wb")

        command = [
            "bgzip",
            "-@",
            str(self.threads),
            "-c",
        ]

        self.process = subprocess.Popen(
            command,
            stdin=subprocess.PIPE,
            stdout=self.output_handle,
            stderr=subprocess.PIPE,
        )

        if self.process.stdin is None:
            raise RuntimeError("Could not open stdin pipe to bgzip.")

        return self

    def write(self, text: str):
        if self.process is None or self.process.stdin is None:
            raise RuntimeError("BGZF writer is not open.")

        try:
            self.process.stdin.write(text.encode("utf-8"))
        except BrokenPipeError:
            stderr = ""

            if self.process.stderr is not None:
                stderr = self.process.stderr.read().decode(
                    "utf-8",
                    errors="replace",
                )

            raise RuntimeError(
                "bgzip terminated unexpectedly.\n"
                f"{stderr}"
            )

    def __exit__(self, exc_type, exc_value, traceback):
        if self.process is None:
            return False

        if self.process.stdin is not None:
            try:
                self.process.stdin.close()
            except BrokenPipeError:
                pass

        stderr_data = b""

        if self.process.stderr is not None:
            stderr_data = self.process.stderr.read()

        return_code = self.process.wait()

        if self.output_handle is not None:
            self.output_handle.close()

        # If Python itself raised an exception, don't replace it with
        # a secondary bgzip error unless necessary.
        if exc_type is not None:
            return False

        if return_code != 0:
            stderr = stderr_data.decode("utf-8", errors="replace")

            raise RuntimeError(
                f"bgzip failed with exit code {return_code}.\n"
                f"{stderr}"
            )

        return False


# =============================================================================
# GT handling
# =============================================================================

def is_het_gt(sample: str) -> bool:
    """
    Determine whether the GT field is heterozygous.

    Examples flagged:
        0/1
        1/0
        0|1
        0/2
        1/2
        1|2

    Examples not flagged:
        0
        1
        0/0
        1/1
        .
        ./.
        .|.

    GT is expected to be the first FORMAT/sample subfield, as required
    by the VCF specification when GT is present.
    """

    if not sample:
        return False

    # Avoid splitting all FORMAT fields.
    colon = sample.find(":")

    if colon == -1:
        gt = sample
    else:
        gt = sample[:colon]

    if "/" in gt:
        alleles = gt.split("/")
    elif "|" in gt:
        alleles = gt.split("|")
    else:
        # Haploid GT or missing genotype.
        return False

    observed = set()

    for allele in alleles:
        if allele != "." and allele != "":
            observed.add(allele)

            # No need to continue once two alleles have been observed.
            if len(observed) > 1:
                return True

    return False


def set_gt_missing(sample: str, haploid_missing: bool = False) -> str:
    """
    Replace GT with missing while keeping every other FORMAT value.

    Default:
        0/1:10,8:18
          ->
        ./.:10,8:18

    With haploid_missing=True:
        0/1:10,8:18
          ->
        .:10,8:18
    """

    missing = "." if haploid_missing else "./."

    colon = sample.find(":")

    if colon == -1:
        return missing

    return missing + sample[colon:]


# =============================================================================
# AD handling
# =============================================================================

def get_ad_index(
    format_string: str,
    cache: Dict[str, Optional[int]],
) -> Optional[int]:
    """
    Return the position of AD within FORMAT.

    FORMAT strings repeat extensively within VCFs, so cache the result.
    """

    cached = cache.get(format_string, "__NOT_CACHED__")

    if cached != "__NOT_CACHED__":
        return cached

    format_fields = format_string.split(":")

    try:
        index = format_fields.index("AD")
    except ValueError:
        index = None

    cache[format_string] = index

    return index


def parse_ad(
    sample: str,
    ad_index: int,
) -> Optional[List[int]]:
    """
    Extract AD values from a sample field.

    Returns:
        [REF, ALT1, ALT2, ...]

    Examples:
        sample = "0/1:23,17:40"
        AD index = 1
        -> [23, 17]

        sample = "0/2:20,3,15:38"
        -> [20, 3, 15]
    """

    # split() is needed here because AD can be located anywhere in FORMAT.
    values = sample.split(":")

    if ad_index >= len(values):
        return None

    raw_ad = values[ad_index]

    if not raw_ad or raw_ad == ".":
        return None

    parts = raw_ad.split(",")

    if len(parts) < 2:
        return None

    depths = []

    for value in parts:
        if value == "" or value == ".":
            return None

        try:
            depth = int(value)
        except ValueError:
            return None

        depths.append(depth)

    return depths


def ad_is_balanced(
    sample: str,
    ad_index: int,
    low: float,
    high: float,
) -> bool:
    """
    Determine whether an ALT allele has an AD depth balanced against REF.

    Rule:
        low <= ALT / REF <= high

    Multi-allelic sites:
        Every ALT is tested against REF individually.

    IMPORTANT:
        If ANY AD value equals zero, the entire sample/position is ignored,
        preserving the requested behaviour.

    Example:
        AD = 20,12
        ratio = 12/20 = 0.60
        -> flagged if low=0.2 and high=1.8

        AD = 20,0,12
        -> NOT flagged because one AD value is zero
    """

    depths = parse_ad(sample, ad_index)

    if depths is None:
        return False

    # Requested zero-depth rule.
    if 0 in depths:
        return False

    ref = depths[0]

    if ref <= 0:
        return False

    for alt in depths[1:]:
        ratio = alt / ref

        if low <= ratio <= high:
            return True

    return False


# =============================================================================
# Input sample list
# =============================================================================

def read_haploid_samples(path: str) -> set:
    """
    Read haploid sample IDs, ignoring blank lines.
    """

    samples = set()

    with open(path, "r") as handle:
        for line in handle:
            sample = line.strip()

            if sample:
                samples.add(sample)

    return samples


# =============================================================================
# --matt mode
# =============================================================================

def write_positions(
    args,
    wanted_samples: set,
):
    """
    Stream through VCF and write flagged CHROM/POS positions.

    No variants are stored in memory.

    As soon as one haploid individual flags a variant, the remaining
    individuals are skipped because --matt only needs one output record
    per genomic position.
    """

    if args.AD:
        output = make_output_filename(
            args.vcf,
            "_AD_positions.txt",
        )
    else:
        output = make_output_filename(
            args.vcf,
            "_het_positions.txt",
        )

    sample_indices = None
    found_samples = []
    ad_cache: Dict[str, Optional[int]] = {}

    variants = 0
    flagged = 0

    with open_maybe_gzip(args.vcf, "rt") as infile, \
            open(output, "w") as outfile:

        outfile.write("#CHROM\tPOS\n")

        for line in infile:

            if line.startswith("##"):
                continue

            if line.startswith("#CHROM"):
                columns = line.rstrip("\r\n").split("\t")

                if len(columns) < 10:
                    raise RuntimeError(
                        "VCF does not contain sample columns."
                    )

                sample_indices = []
                found_samples = []

                for index in range(9, len(columns)):
                    sample = columns[index]

                    if sample in wanted_samples:
                        sample_indices.append(index)
                        found_samples.append(sample)

                if not sample_indices:
                    raise RuntimeError(
                        "None of the individuals in the haploid list "
                        "were found in the VCF."
                    )

                print(
                    f"Found {len(found_samples)} haploid individual(s) "
                    f"in the VCF.",
                    file=sys.stderr,
                )

                continue

            if line.startswith("#"):
                continue

            if sample_indices is None:
                raise RuntimeError(
                    "Could not find the #CHROM VCF header line."
                )

            variants += 1

            fields = line.rstrip("\r\n").split("\t")

            # Protect against malformed rows.
            if len(fields) <= max(sample_indices):
                continue

            is_flagged = False

            if args.AD:
                if len(fields) < 9:
                    continue

                format_string = fields[8]

                ad_index = get_ad_index(
                    format_string,
                    ad_cache,
                )

                if ad_index is not None:
                    for column_index in sample_indices:
                        if ad_is_balanced(
                            fields[column_index],
                            ad_index,
                            args.low,
                            args.high,
                        ):
                            is_flagged = True

                            # Major optimization for --matt:
                            # one positive individual is sufficient.
                            break

            else:
                for column_index in sample_indices:
                    if is_het_gt(fields[column_index]):
                        is_flagged = True
                        break

            if is_flagged:
                outfile.write(
                    f"{fields[0]}\t{fields[1]}\n"
                )
                flagged += 1

            if (
                args.progress > 0
                and variants % args.progress == 0
            ):
                print(
                    f"Processed {variants:,} variants; "
                    f"flagged {flagged:,}",
                    file=sys.stderr,
                )

    print(
        f"Processed {variants:,} variants.",
        file=sys.stderr,
    )
    print(
        f"Flagged {flagged:,} positions.",
        file=sys.stderr,
    )
    print(
        f"Positions saved as: {output}",
        file=sys.stderr,
    )


# =============================================================================
# Modified VCF mode
# =============================================================================

def write_modified_vcf(
    args,
    wanted_samples: set,
):
    """
    Stream the VCF, modify selected sample GTs, and write proper
    BGZF-compressed VCF output.
    """

    output = make_output_filename(
        args.vcf,
        "_modified.vcf.gz",
    )

    sample_indices = None
    found_samples = []

    ad_cache: Dict[str, Optional[int]] = {}

    variants = 0
    modified_genotypes = 0
    modified_positions = 0

    with open_maybe_gzip(args.vcf, "rt") as infile, \
            BgzipWriter(output, threads=args.threads) as outfile:

        for line in infile:

            # -------------------------------------------------------------
            # VCF metadata
            # -------------------------------------------------------------

            if line.startswith("##"):
                outfile.write(line)
                continue

            # -------------------------------------------------------------
            # Main VCF header
            # -------------------------------------------------------------

            if line.startswith("#CHROM"):
                columns = line.rstrip("\r\n").split("\t")

                if len(columns) < 10:
                    raise RuntimeError(
                        "VCF does not contain sample columns."
                    )

                sample_indices = []
                found_samples = []

                for index in range(9, len(columns)):
                    sample = columns[index]

                    if sample in wanted_samples:
                        sample_indices.append(index)
                        found_samples.append(sample)

                if not sample_indices:
                    raise RuntimeError(
                        "None of the individuals in the haploid list "
                        "were found in the VCF."
                    )

                print(
                    f"Found {len(found_samples)} haploid individual(s) "
                    f"in the VCF.",
                    file=sys.stderr,
                )

                outfile.write(line)

                continue

            # Preserve any unusual additional header lines.
            if line.startswith("#"):
                outfile.write(line)
                continue

            if sample_indices is None:
                raise RuntimeError(
                    "Could not find the #CHROM VCF header line."
                )

            # -------------------------------------------------------------
            # Variant
            # -------------------------------------------------------------

            variants += 1

            fields = line.rstrip("\r\n").split("\t")

            if len(fields) <= max(sample_indices):
                # Preserve malformed/unexpected rows rather than silently
                # discarding them.
                outfile.write(line)
                continue

            position_modified = False

            # -------------------------------------------------------------
            # AD MODE
            # -------------------------------------------------------------

            if args.AD:
                if len(fields) >= 9:
                    format_string = fields[8]

                    ad_index = get_ad_index(
                        format_string,
                        ad_cache,
                    )

                    if ad_index is not None:

                        for column_index in sample_indices:

                            sample = fields[column_index]

                            if ad_is_balanced(
                                sample,
                                ad_index,
                                args.low,
                                args.high,
                            ):
                                fields[column_index] = set_gt_missing(
                                    sample,
                                    haploid_missing=args.haploid_missing,
                                )

                                modified_genotypes += 1
                                position_modified = True

            # -------------------------------------------------------------
            # GT MODE
            # -------------------------------------------------------------

            else:

                for column_index in sample_indices:

                    sample = fields[column_index]

                    if is_het_gt(sample):
                        fields[column_index] = set_gt_missing(
                            sample,
                            haploid_missing=args.haploid_missing,
                        )

                        modified_genotypes += 1
                        position_modified = True

            if position_modified:
                modified_positions += 1

            outfile.write(
                "\t".join(fields) + "\n"
            )

            # -------------------------------------------------------------
            # Lightweight progress reporting
            # -------------------------------------------------------------

            if (
                args.progress > 0
                and variants % args.progress == 0
            ):
                print(
                    f"Processed {variants:,} variants; "
                    f"modified {modified_genotypes:,} genotypes "
                    f"at {modified_positions:,} positions.",
                    file=sys.stderr,
                )

    print(
        f"Processed {variants:,} variants.",
        file=sys.stderr,
    )

    print(
        f"Modified {modified_genotypes:,} genotype(s) "
        f"at {modified_positions:,} position(s).",
        file=sys.stderr,
    )

    print(
        f"BGZF-compressed VCF saved as: {output}",
        file=sys.stderr,
    )


# =============================================================================
# Main
# =============================================================================

def main():

    parser = argparse.ArgumentParser(
        description=(
            "Set heterozygous/balanced positions of specified haploid "
            "individuals to missing in a VCF. Uses streaming processing "
            "and writes proper BGZF-compressed VCF output."
        )
    )

    parser.add_argument(
        "-v",
        "--vcf",
        required=True,
        help=(
            "Input VCF. Supported: .vcf, .vcf.gz, .vcf.bgz, .vcf.bgzip"
        ),
    )

    parser.add_argument(
        "-l",
        "--list",
        required=True,
        help=(
            "Text file containing haploid sample IDs, one per line."
        ),
    )

    # Retained so old command lines using -r do not fail.
    # It is deliberately ignored because header detection is automatic.
    parser.add_argument(
        "-r",
        "--rownum",
        type=int,
        required=False,
        help=(
            "Deprecated and ignored. VCF header rows are now detected "
            "automatically."
        ),
    )

    parser.add_argument(
        "--matt",
        action="store_true",
        help=(
            "Only output flagged CHROM/POS positions. "
            "Do not create a modified VCF."
        ),
    )

    parser.add_argument(
        "--AD",
        action="store_true",
        help=(
            "Use AD allele-depth balance instead of GT heterozygosity."
        ),
    )

    parser.add_argument(
        "--low",
        type=float,
        default=0.2,
        help=(
            "Lower inclusive ALT/REF AD ratio in --AD mode "
            "(default: 0.2)."
        ),
    )

    parser.add_argument(
        "--high",
        type=float,
        default=1.8,
        help=(
            "Upper inclusive ALT/REF AD ratio in --AD mode "
            "(default: 1.8)."
        ),
    )

    parser.add_argument(
        "--threads",
        type=int,
        default=4,
        help=(
            "Number of bgzip compression threads "
            "(default: 4)."
        ),
    )

    parser.add_argument(
        "--progress",
        type=int,
        default=1_000_000,
        help=(
            "Print progress every N variants. "
            "Use 0 to disable progress messages "
            "(default: 1000000)."
        ),
    )

    parser.add_argument(
        "--haploid-missing",
        action="store_true",
        help=(
            "Write missing GT as '.' instead of './.'. "
            "Use this for strictly haploid VCF representation."
        ),
    )

    args = parser.parse_args()

    # ---------------------------------------------------------------------
    # Validate options
    # ---------------------------------------------------------------------

    if args.low < 0:
        parser.error("--low cannot be negative.")

    if args.high < args.low:
        parser.error("--high must be greater than or equal to --low.")

    if args.threads < 1:
        parser.error("--threads must be at least 1.")

    if args.progress < 0:
        parser.error("--progress cannot be negative.")

    if not os.path.isfile(args.vcf):
        parser.error(
            f"VCF file does not exist: {args.vcf}"
        )

    if not os.path.isfile(args.list):
        parser.error(
            f"Haploid sample list does not exist: {args.list}"
        )

    if not args.matt and shutil.which("bgzip") is None:
        parser.error(
            "bgzip was not found in PATH. "
            "Install HTSlib/tabix/bcftools before writing BGZF output."
        )

    # ---------------------------------------------------------------------
    # Read haploid sample IDs
    # ---------------------------------------------------------------------

    wanted_samples = read_haploid_samples(args.list)

    if not wanted_samples:
        parser.error(
            "The haploid individual list is empty."
        )

    print(
        f"Loaded {len(wanted_samples)} sample ID(s) "
        f"from {args.list}.",
        file=sys.stderr,
    )

    # ---------------------------------------------------------------------
    # Process
    # ---------------------------------------------------------------------

    if args.matt:
        write_positions(
            args,
            wanted_samples,
        )
    else:
        write_modified_vcf(
            args,
            wanted_samples,
        )


if __name__ == "__main__":
    try:
        main()

    except KeyboardInterrupt:
        print(
            "\nInterrupted by user.",
            file=sys.stderr,
        )
        sys.exit(130)

    except Exception as exc:
        print(
            f"\nERROR: {exc}",
            file=sys.stderr,
        )
        sys.exit(1)
