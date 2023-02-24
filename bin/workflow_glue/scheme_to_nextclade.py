#!/usr/bin/env python
"""Script to convert ARTIC primer scheme bed to nextclade input."""

import os

import pysam

from .util import get_named_logger, wf_parser  # noqa: ABS101


def main(args):
    """Run entry point."""
    logger = get_named_logger("nextclade-prep")

    pysam.FastaFile(args.reference)
    refs = {r.name: r.sequence for r in pysam.FastxFile(args.reference)}
    scheme = os.path.splitext(os.path.basename(args.bed))[0]
    # MN908947.3    30    54    nCoV-2019_1_LEFT    1    +
    with open(args.bed) as fh, open(args.output, 'w') as out_fh:
        out_fh.write("Country (Institute),Target,Oligonucleotide,Sequence\n")
        for line in fh.readlines():
            # the V1200 bed doesn't specific orientation
            rname, start, end, name, pool, *_ = line.split('\t')
            pname = '{}:{}-{}'.format(rname, start, end)
            seq = refs[rname][int(start):int(end)]
            out_fh.write(','.join((scheme, pname, name, seq)))
            out_fh.write('\n')

    logger.info(f"Nextclade scheme written to {args.output}.")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("nextclade-prep")

    parser.add_argument('bed')
    parser.add_argument('reference')
    parser.add_argument('output')

    return parser
