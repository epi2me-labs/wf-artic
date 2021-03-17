#!/usr/bin/env python
"""Script to generate simple per read metrics."""

import argparse
import glob
import itertools
import os

import numpy as np
import pysam


def mean_qual(quals):
    """Calculate mean quality score of a read."""
    qual = np.fromiter(
        (ord(x) - 33 for x in quals),
        dtype=int, count=len(quals))
    mean_p = np.mean(np.power(10, qual / -10))
    return -10 * np.log10(mean_p)


def main():
    """Run entry point."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "directory", help="Directory containing .fastq(.gz) files")
    parser.add_argument(
        "output", help="Output file")
    parser.add_argument(
        "--sample_name", help="Sample name to add as column to output")
    args = parser.parse_args()

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


if __name__ == "__main__":
    main()
