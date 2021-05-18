#!/usr/bin/env python
"""A simple left-join genotyper."""

import argparse
from datetime import datetime

import pandas as pd
import pysam

CHROM = 'CHROM'
POS = 'POS'
ID = 'ID'
REF = 'REF'
ALT = 'ALT'
QUAL = 'QUAL'
FILTER = 'FILTER'
INFO = 'INFO'

CALLS = 'CALLS'
COV = 'COVERAGE'

vcf_cols = [CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO]


def read_vcf(fname):
    """Read a VCF file."""
    df = pd.read_csv(fname, comment='#', header=None, delimiter='\t')
    df = df.rename(columns=dict(enumerate(vcf_cols)))
    return df


def extract_info(record):
    """Expand the info field for the given variant call record."""
    extracted = {}
    for item in record.split(';'):
        item_split = item.split('=')
        extracted[item_split[0]] = item_split[1]

    return extracted


def determine_coverage(bam, chrom, pos):
    """Get read coverage at the given coordinates."""
    try:
        column = next(bam.pileup(
            chrom, pos, pos + 1,
            truncate=True)
        )
    except Exception:
        coverage = 0
    else:
        coverage = column.n
    return coverage


def extract_calls(df, chrom, pos, ref):
    """Get variant calls for the given coordinates."""
    df_calls = df[
        (df[CHROM] == chrom)
        & (df[POS] == pos)
        & (df[REF] == ref)
    ]

    if df_calls.empty:
        return []

    df_calls[QUAL] = pd.to_numeric(df_calls[QUAL])

    calls_list = []
    for call in df_calls.to_dict('records'):
        calls_list.append({
            ALT: call[ALT],
            QUAL: call[QUAL],
            INFO: extract_info(call[INFO])
        })

    # Retain the best call by QUAL for this locus
    calls = sorted(calls_list, key=lambda k: k[QUAL])
    return calls


def get_results(vcf_file, variants_file, bam_file):
    """Acquire the results of variant calling wrt the variants file."""
    results = []
    bam = pysam.AlignmentFile(bam_file, "rb")

    # Load in data or an empty frame
    try:
        df_vcf = read_vcf(vcf_file)
    except Exception:
        df_vcf = pd.DataFrame(columns=vcf_cols)

    # Read varants file
    df_variants = read_vcf(variants_file)
    for _, row in df_variants.iterrows():
        data = {
            CHROM: row[CHROM],
            POS: row[POS],
            ID: row[ID],
            REF: row[REF],
            ALT: row[ALT],
            CALLS: extract_calls(
                df_vcf, row[CHROM],
                row[POS], row[REF]
            ),
            COV: determine_coverage(
                bam, row[CHROM], row[POS]
            )
        }

        results.append(data)
    return results


def compile_table(results, sample, lab_id, timestamp, testkit,
                  valid_coverage=20, valid_quality=20):
    """Compile the .csv output from the variant calling results."""
    table_data = []
    for entry in results:
        # Begin with default values
        result = "No mutation"
        CH1_result = "wt"
        CH1_target = entry[ID]
        CH1_conf = 0

        # If a call exists for this locus
        if calls := entry[CALLS]:
            top_call = calls[0]
            # Note this is a simplification of the likelihood
            # assigned to this call because it does not take
            # into account the probabilities given to any
            # of the alternatives.
            CH1_conf = 1 - (10 ** (-top_call[QUAL] / 10))
            CH1_conf = f"{CH1_conf:.2f}"
            # If the call is to the expected alt
            # we consider it an expected "mutation"
            # Note, we do not expect bi-allelic calls
            # to appear in medaka output
            if top_call['ALT'] == entry['ALT']:
                CH1_result = "mt"
                result = "Mutation"
            # Else if the call is not the ref (and also not
            # the expected alt), then we are provisionally
            # marking the call as "variant"
            elif top_call['ALT'] not in ['.', entry['REF']]:
                CH1_result = "vr"
                result = "Variant"
            # Independent of the above, if the quality of
            # the call does not meet the minimum threshold,
            # mark the call as invalid
            if top_call[QUAL] < valid_quality:
                CH1_result = "n/a"
                result = "Low quality"

        # Filter by coverage
        if entry[COV] < valid_coverage:
            CH1_result = "n/a"
            result = "No Amplification"
            CH1_conf = 0

        table_data.append({
            'Sample': sample,
            'Result': result,
            'Date Tested': timestamp,
            'Lab ID': lab_id,
            'testKit': testkit,
            'CH1-Target': CH1_target,
            'CH1-Result': CH1_result,
            'CH1-Conf': CH1_conf
        })

    table_df = pd.DataFrame(table_data)
    table_df = table_df.sort_values(
        by=['Sample', 'CH1-Target'], ascending=False)
    return table_df


def main():
    """Create the genotype summary output the result to file."""
    parser = argparse.ArgumentParser()
    tz = datetime.utcnow().astimezone().tzinfo
    date = datetime.today().strftime('%Y-%m-%d %I:%M:%S %p')
    timestamp = date + f" {tz}"

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
        help="Path to outfile")
    parser.add_argument(
        "--valid_coverage", type=int, default=20,
        help="Coverage required for valid call.")
    # Meta fields
    parser.add_argument(
        "-s", "--sample", required=False, type=str,
        help="Override the default sample name")
    parser.add_argument(
        "--timestamp", required=False, type=str,
        help="Override the automatically generated timestamp",
        default=timestamp)
    parser.add_argument(
        "--lab_id", required=False, type=str,
        help="Override the default lab_id (ONT)", default="ONT")
    parser.add_argument(
        "--testkit", required=False, type=str,
        help="Override default testKit (ONT_spikeseq)",
        default="ONT_spikeseq")

    args = parser.parse_args()

    sample = args.sample
    if not sample:
        sample = args.vcf.split('.')[0]

    results = get_results(args.vcf, args.variants, args.bam)
    table = compile_table(results, sample, args.lab_id,
                          args.timestamp, args.testkit)

    table.to_csv(args.outfile, index=False, na_rep='n/a')


if __name__ == "__main__":
    main()
