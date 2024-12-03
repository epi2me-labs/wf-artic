#!/usr/bin/env python
"""Create report file."""

from aplanat import report
from aplanat.components import simple as scomponents

from .util import get_named_logger, wf_parser  # noqa: ABS101


def main(args):
    """Run entrypoint."""
    logger = get_named_logger("report-error")

    report_doc = report.WFReport(
        "SARS-CoV-2 ARTIC Sequencing report", "wf-artic",
        revision=args.revision, commit=args.commit)
    section = report_doc.add_section()
    # error message
    section.markdown('###**' + args.error_message + '**')
    # Versions and params
    section = report_doc.add_section(
        section=scomponents.version_table(args.versions))
    section = report_doc.add_section(
        section=scomponents.params_table(args.params))

    # write report
    report_doc.write(args.output)

    logger.info(f"Written report to {args.output}.")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("report-error")

    parser.add_argument(
        "--output",
        help="Report output filename")
    parser.add_argument(
        "--revision", default='unknown',
        help="git branch/tag of the executed workflow")
    parser.add_argument(
        "--commit", default='unknown',
        help="git commit of the executed workflow")
    parser.add_argument(
        "--params", default=None,
        help="A csv containing the parameter key/values")
    parser.add_argument(
        "--versions",
        help="directory contained CSVs containing name,version.")
    parser.add_argument(
        "--error_message",
        help="directory contained CSVs containing name,version.")

    return parser
