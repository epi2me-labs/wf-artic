#!/usr/bin/env python
"""A simple left-join genotyper."""

import argparse

import pandas as pd
import pysam


vcf_cols = "CHROM POS ID REF ALT QUAL FILTER INFO".split()


def read_vcf(fname):
    """Read a VCF file."""
    df = pd.read_csv(fname, comment='#', header=None, delimiter='\t')
    df = df.rename(columns=dict(enumerate(vcf_cols)))
    return df


def create_genotype_summary(vcf_file, bam_file, variants_file, valid_coverage):
    """Merge the variants tsv file with the calls in the vcf."""
    # Grab sample name from vcf filename
    sample = vcf_file.split('.')[0]  # TODO: this is a bit flakey

    # Convert the calls from artic to a vcf
    try:
        df = read_vcf(vcf_file)
    except Exception:
        df = pd.DataFrame(columns=vcf_cols)
    df['variant_score'] = df['QUAL']

    # Read in the variants we have stored as a tsv
    # and add some meta fields
    df_variants = read_vcf(variants_file)
    df_variants['variant'] = df_variants['ID']

    # Merge the variants file with the vcf from artic
    # retaining the variants file overlaps
    df_merged = df_variants.merge(
        df, how='left', on=['CHROM', 'POS', 'REF', 'ALT'])
    # the ARTIC VCF can have variants multiple times (from multiple PCR pools)
    df_merged = df_merged \
        .sort_values('variant_score', ascending=False) \
        .drop_duplicates('variant')

    # Re-write sample name
    df_merged['sample'] = sample

    # Adjust fields depending on whether there is
    # a mutation present or not
    df_merged['result'] = 'Not present'
    df_merged['coverage'] = 0
    df_merged['status'] = 'Valid'
    df_merged.loc[~df_merged['variant_score'].isna(), 'result'] = 'Present'
    df_merged.loc[df_merged['variant_score'].isna(), 'variant_score'] = 'N/A'

    # Amend status column based on coverage
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for index, row in df_merged.iterrows():
            try:
                column = next(bam.pileup(
                    row["CHROM"], row["POS"], row["POS"] + 1, truncate=True))
            except Exception:
                coverage = 0
            else:
                coverage = column.n
            df_merged.loc[index, 'coverage'] = coverage
            if coverage < valid_coverage:
                status = "No data" if coverage == 0 else "Low coverage"
                df_merged.loc[index, 'status'] = status

    return df_merged[
        ['sample', 'variant', 'variant_score', 'coverage', 'result', 'status']
    ].sort_values(['sample', 'variant'])


def main():
    """Create the genotype summary output the result to file."""

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-v", "--vcf", required=True, type=str,
        help="Path to VCF file")
    parser.add_argument(
        "-b", "--bam", required=True, type=str,
        help="Path to BAM file")
    parser.add_argument(
        "-d", "--variants", required=True, type=str,
        help="Path to variants of interest file")
    parser.add_argument(
        "-o", "--outfile",
        help="Path  to outfile")
    parser.add_argument(
        "--valid_coverage", type=int, default=20,
        help="Coverage required for valid call.")

    args = parser.parse_args()

    result = create_genotype_summary(
        args.vcf, args.bam, args.variants, args.valid_coverage)
    result.to_csv(args.outfile, index=False)


if __name__ == "__main__":
    main()
