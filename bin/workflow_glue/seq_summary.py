#!/usr/bin/env python
"""Script to generate simple per read metrics."""

import glob
import itertools
import os

import numpy as np
import pysam

from .util import get_named_logger, wf_parser  # noqa: ABS101


def mean_qual(quals):
    """Calculate mean quality score of a read."""
    qual = np.fromiter(
        (ord(x) - 33 for x in quals),
        dtype=int, count=len(quals))
    mean_p = np.mean(np.power(10, qual / -10))
    return -10 * np.log10(mean_p)


def main(args):
    """Run entrypoint."""
    logger = get_named_logger("seq-summary")

    fastqs = glob.glob(os.path.join(args.directory, "*.fastq*"))
    reads = itertools.chain.from_iterable(
        pysam.FastxFile(fname) for fname in fastqs)

    with open(args.output, "w") as fh:
        # names as in Guppy
        fh.write(
            "sample_name\tread_id\tsequence_length_template\t"
            "mean_qscore_template\n")
        for read in reads:
            fh.write(
                "\t".join(
                    str(x) for x in (
                        args.sample_name, read.name, len(read.sequence),
                        mean_qual(read.quality))))
            fh.write("\n")

    logger.info(f"Seq summary written to {args.output}.")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("seq-summary")

    parser.add_argument(
        "directory", help="Directory containing .fastq(.gz) files")
    parser.add_argument(
        "output", help="Output file")
    parser.add_argument(
        "--sample_name", help="Sample name to add as column to output")

    return parser
