#!/usr/bin/env python
"""Create telemetry file."""

import argparse
import json
import statistics
import sys

import pandas as pd
from pysam import AlignmentFile, FastxFile, VariantFile


# Reporting
VERSION = 'version'
REFERENCE = 'reference'
BARCODES = 'barcodes'
AMPLICONS = 'amplicons'
PRIMERS = 'primers'
COVERAGE = 'coverage'
CALLS = 'calls'
ACTUAL = 'actual'
OVERLAP = 'overlap'

# BED
START = 'start'
END = 'end'
PRIMER = 'primer'
PRIMER_SET = 'primer_set'
STRAND = 'strand'
POOL = 'pool'

# VCF
CHROM = 'chrom'
POS = 'pos'
QUAL = 'qual'
REF = 'ref'
ALT = 'alt'

# Stats
MEAN = 'mean'
MEDIAN = 'median'
READS = 'reads'
ABOVE_20X = 'perc_bases_above_20x'
ABOVE_150X = 'perc_bases_above_150x'


def get_reference_name(reference, index):
    """Get name of first sequence in fastx file."""
    with FastxFile(reference) as ref:
        for idx, entry in enumerate(ref):
            if idx != index:
                continue
            return entry.name


def get_scheme(scheme_bedfile):
    """Load primer scheme bedfile."""
    return pd.read_csv(
        scheme_bedfile,
        sep="\t",
        names=[CHROM, START, END, PRIMER, PRIMER_SET, STRAND])


def groupby_amplicons(df, ind):
    """Split primer name to retain amplicon name."""
    splits = df[PRIMER].loc[ind].split('_')
    return '_'.join([splits[0], splits[1]])


def get_alignment_tag(alignment, tag):
    """Get SAM tag."""
    return alignment.get_tag(tag)


def get_median(arr, default=0):
    """Get median value of array."""
    try:
        return "{:.2f}".format(statistics.median(arr))
    except statistics.StatisticsError:
        return default


def get_mean(arr, default=0):
    """Get mean value of array."""
    try:
        return "{:.2f}".format(statistics.mean(arr))
    except statistics.StatisticsError:
        return default


def build_telemetry(
        scheme_name, scheme_bed,
        reference, samples, alignments, calls):
    """Return telemetry blob."""
    reference_name = get_reference_name(reference, 0)

    telemetry = {
        VERSION: {
            PRIMER_SET: scheme_name,
            REFERENCE: reference_name,
        },
        BARCODES: {}
    }

    # Minor check that all sample data is present
    for arg in [alignments, calls]:
        if len(arg) != len(samples):
            print(
                "Error: samples, alignments and calls "
                "must be the same length.")
            sys.exit(1)

    # Load scheme data and get amplicon groups
    df_scheme = get_scheme(scheme_bed)
    grouped_by_amplicon = df_scheme.groupby(
        lambda x: groupby_amplicons(df_scheme, x))

    # Iterate over samples and construct data
    sample_data = list(zip(samples, alignments, calls))
    for sample_data in sample_data:
        telemetry[BARCODES][sample_data[0]] = {}
        sample = telemetry[BARCODES][sample_data[0]]

        bam = AlignmentFile(sample_data[1])
        vcf = VariantFile(sample_data[2])

        sample[AMPLICONS] = {}
        for amp_name, primers in grouped_by_amplicon:
            primers.reset_index(inplace=True)
            sample[AMPLICONS][amp_name] = {}
            amplicon = sample[AMPLICONS][amp_name]

            # Gather amplicon level data
            pool = primers.at[0, PRIMER_SET]
            start = primers[START].min()
            end = primers[END].max()
            length = end - start

            # Initialise variables to accumulate
            read_ids = set()
            reads_per_base = []
            stats_per_primer = {
                p[0]: {
                    START: p[1], END: p[2], ACTUAL: [],
                    OVERLAP: [], CALLS: []
                }
                for p in primers[[PRIMER, START, END]].values
            }
            bases_above_20x = 0
            bases_above_150x = 0

            # Assuming that coverage has not been tampered with upstream
            for base in bam.pileup(
                    reference_name, start, end, min_base_quality=7):
                # Pileup will return all bases from overlapping reads
                # Ignore those outside the amplicon bounds
                if not start < base.pos < end:
                    continue

                # Count bases and reads at this position
                # For the entire amplicon
                per_base_count = 0
                per_base_count_alt = 0
                for read in base.pileups:
                    rg = get_alignment_tag(read.alignment, 'RG')
                    if int(rg) == int(pool):
                        read_ids.add(read.alignment.query_name)
                        per_base_count += 1
                        continue
                    per_base_count_alt += 1
                reads_per_base.append(per_base_count)

                # And as required for any primer site we're currently within
                for primer in stats_per_primer.values():
                    if primer[START] <= base.pos <= primer[END]:
                        primer[ACTUAL].append(per_base_count)
                        primer[OVERLAP].append(per_base_count_alt)

                        # Cheekily grab calls from overlaps whilst we're at it
                        try:
                            for call in vcf.fetch(
                                    reference_name, base.pos-1, base.pos):
                                primer[CALLS].append({
                                    POS: call.pos,
                                    QUAL: call.qual,
                                    REF: call.ref,
                                    ALT: call.alts,
                                    POOL: call.info.get('Pool')
                                })
                        except ValueError:
                            pass

                # Calculate additional summary stats
                if per_base_count > 20:
                    bases_above_20x += 1
                if per_base_count > 150:
                    bases_above_150x += 1

            # Put all the data together per amplicon
            amplicon[POOL] = int(pool)

            amplicon[COVERAGE] = {
                READS: len(read_ids),
                MEDIAN: get_median(reads_per_base),
                MEAN: get_mean(reads_per_base),
                ABOVE_20X: "{:.2f}".format(bases_above_20x / length * 100),
                ABOVE_150X: "{:.2f}".format(bases_above_150x / length * 100)
            }

            amplicon[PRIMERS] = {
                k: {
                    COVERAGE: {
                        ACTUAL: {
                            MEAN: get_mean(v[ACTUAL]),
                            MEDIAN: get_median(v[ACTUAL]),
                        },
                        OVERLAP: {
                            MEAN: get_mean(v[OVERLAP]),
                            MEDIAN: get_median(v[OVERLAP]),
                        },
                    },
                    CALLS: v[CALLS]
                } for k, v in stats_per_primer.items()
            }

    return telemetry


def main():
    """Run entry point."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "output",
        help="Telemetry output filename")
    parser.add_argument(
        "--scheme_name", required=True,
        help="Primer schema name")
    parser.add_argument(
        "--scheme_bed", required=True,
        help="Primer schema BED file")
    parser.add_argument(
        "--reference", required=True,
        help="Reference file used with artic")
    parser.add_argument(
        "--samples", nargs='+', required=True,
        help="Sample names")
    parser.add_argument(
        "--calls", nargs='+', required=True,
        help="Variant calls made by artic")
    parser.add_argument(
        "--alignments", nargs='+', required=True,
        help="Alignments made by artic")
    args = parser.parse_args()

    telemetry = build_telemetry(
        scheme_name=args.scheme_name,
        scheme_bed=args.scheme_bed,
        reference=args.reference,
        samples=args.samples,
        alignments=args.alignments,
        calls=args.calls)

    # Write to JSON
    with open("telemetry.json", 'w') as fout:
        data = json.dumps(telemetry, indent=4, separators=(',', ':'))
        fout.write(data)


if __name__ == "__main__":
    main()
