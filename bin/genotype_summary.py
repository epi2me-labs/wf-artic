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


def create_genotype_summary(vcf_file, bam_file, variants_file):
    """Merge the variants tsv file with the calls in the vcf."""
    # Grab sample name from vcf filename
    sample = vcf_file.split('.')[0]

    # Convert the calls from artic to a vcf
    try:
        df = read_vcf(vcf_file)
        df['variant_score'] = df['QUAL']
    except Exception:
        df = pd.DataFrame(columns=vcf_cols)

    # Read in the variants we have stored as a tsv
    # and add some meta fields
    df_variants = read_vcf(variants_file)
    df_variants['variant'] = df_variants['ID']

    # Merge the variants file with the vcf from artic
    # retaining the variants file overlaps
    df_merged = df_variants.merge(
        df, how='left', on=['CHROM', 'POS', 'REF', 'ALT'])

    # Re-write sample name
    df_merged['sample'] = sample

    # Adjust fields depending on whether there is
    # a mutation present or not
    df_merged['result'] = 'no_mutation'
    df_merged.loc[~df_merged['variant_score'].isna(), 'result'] = 'mutation'
    df_merged.loc[df_merged['variant_score'].isna(), 'variant_score'] = 'N/A'

    # Iterate over the merged variants and attempt to
    # find the coverage for each retained variant
    # TODO: this could go, the VCF from artic contains already depth info
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for index, row in df_merged.iterrows():
            try:
                covered = len(list(bam.pileup(
                    row['CHROM'], row["POS"], row["POS"] + 1, truncate=True)))
            except ValueError:
                covered = None
            if not covered:
                df_merged.loc[index, 'result'] = "no_amplification"

    return df_merged[
        ['sample', 'variant', 'variant_score', 'result']
    ].sort_values(['sample', 'variant'])


def main(vcf_file, bam_file, variants_file, outfile):
    """Create the genotype summary output the result to file."""
    result = create_genotype_summary(
        vcf_file,
        bam_file,
        variants_file,
    )

    result.to_csv(outfile, index=False)


def parse_args():
    """Validate command line arguments."""
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-v", "--vcf", required=True, type=str,
        help="Path to 50x coverage filtered vcf"
    )
    parser.add_argument(
        "-b", "--bam", required=True, type=str,
        help="Path to bam file"
    )
    parser.add_argument(
        "-d", "--variants", required=True, type=str,
        help="Path to variants of interest file"
    )
    parser.add_argument(
        "-o", "--outfile", required=True, type=str,
        help="Path  to outfile"
    )

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    main(
        args.vcf,
        args.bam,
        args.variants,
        args.outfile
    )
