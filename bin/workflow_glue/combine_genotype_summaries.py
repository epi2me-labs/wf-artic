#!/usr/bin/env python
"""Create report file."""

import pandas as pd

from .util import get_named_logger, wf_parser  # noqa: ABS101


def read_files(summaries, sep='\t'):
    """Read a set of files and join to single dataframe."""
    dfs = list()
    for fname in sorted(summaries):
        dfs.append(pd.read_csv(fname, sep=sep))
    return pd.concat(dfs)


def main(args):
    """Run entrypoint."""
    logger = get_named_logger("combine-genotype-summaries")

    df = read_files(args.genotypes, sep=',')
    df = df[[
        'Sample', 'Result', 'Date Tested', 'Lab ID', 'testKit',
        'CH1-Target', 'CH1-Result', 'CH1-Conf'
    ]]
    df = df.sort_values(
        by=['Sample', 'CH1-Target'],
        ascending=True
    )

    df.to_csv(args.output, index=False, na_rep='n/a')

    logger.info(f"Written genotype summaries to {args.output}.")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("combine-genotype-summaries")
    parser.add_argument(
        "-o", "--output", required=True,
        help="Output filename")
    parser.add_argument(
        "-g", "--genotypes", nargs='+', required=True,
        help="Genotyping summary files")
    return parser
