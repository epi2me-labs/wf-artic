#!/usr/bin/env python

import argparse
import glob
import os

import aplanat
from aplanat import annot, bars, gridplot, hist, lines, points, report
from aplanat.util import Colors
from aplanat.components import bcfstats
from bokeh.layouts import gridplot, layout
from bokeh.models import Panel, Tabs
import numpy as np
import pandas as pd

def read_files(summaries):
    dfs = list()
    for fname in sorted(summaries):
        dfs.append(pd.read_csv(fname, sep="\t"))
    return pd.concat(dfs)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("nextclade", help="nexclade json output file")
    parser.add_argument("output", help="Report output filename")
    parser.add_argument("--depths", nargs='+', required=True, help="Depth summary files")
    parser.add_argument("--summaries", nargs='+', required=True, help="Sequencing summary files")
    parser.add_argument("--bcftools_stats", nargs='+', required=True, help="Outputs from bcftools stats")
    parser.add_argument("--min_len", default=300, type=int, help="Minimum read length")
    parser.add_argument("--max_len", default=700, type=int, help="Maximum read length")
    args = parser.parse_args()

    report_doc = report.HTMLReport(
        "SARS-CoV-2 ARTIC Sequencing report",
        "Results generated through the wf-artic nextflow workflow provided by Oxford Nanopore Technologies")
    section = report_doc.add_section()

    section.markdown('''
### Read Quality control
This section displays basic QC metrics indicating read data quality.
''')

    # read length summary
    seq_summary = read_files(args.summaries)
    total_bases = seq_summary['sequence_length_template'].sum()
    mean_length = total_bases / len(seq_summary)
    median_length = np.median(seq_summary['sequence_length_template'])
    datas = [seq_summary['sequence_length_template']]
    length_hist = hist.histogram(
        datas, colors=[Colors.cerulean], bins=100,
        title="Read length distribution.",
        x_axis_label='Read Length / bases',
        y_axis_label='Number of reads',
        xlim=(0, 2000))
    length_hist = annot.marker_vline(
        length_hist, args.min_len,
        label="Min: {}".format(args.min_len), text_baseline='bottom', color='grey')
    length_hist = annot.marker_vline(
        length_hist, args.max_len,
        label="Max: {}".format(args.max_len), text_baseline='top')
    length_hist = annot.subtitle(
        length_hist,
        "Mean: {:.0f}. Median: {:.0f}".format(
            mean_length, median_length))

    datas = [seq_summary['mean_qscore_template']]
    mean_q, median_q = np.mean(datas[0]), np.median(datas[0])
    q_hist = hist.histogram(
        datas, colors=[Colors.cerulean], bins=100,
        title="Read quality score",
        x_axis_label="Quality score",
        y_axis_label="Number of reads",
        xlim=(4, 25))
    q_hist = annot.subtitle(
        q_hist,
        "Mean: {:.0f}. Median: {:.0f}".format(
            mean_q, median_q))

    # barcode count plot
    good_reads = seq_summary.loc[
        (seq_summary['sequence_length_template'] > args.min_len)
        & (seq_summary['sequence_length_template'] < args.max_len)]
    barcode_counts = (
        pd.DataFrame(good_reads['sample_name'].value_counts())
        .sort_index()
        .reset_index()
        .rename(
            columns={'index':'sample', 'sample_name':'count'})
        )

    bc_counts = bars.simple_bar(
        barcode_counts['sample'].astype(str), barcode_counts['count'], colors=[Colors.cerulean]*len(barcode_counts),
        title='Number of reads per barcode (filtered by {} < length < {})'.format(args.min_len, args.max_len),
        plot_width=None
    )
    bc_counts.xaxis.major_label_orientation = 3.14/2
    section.plot(
        layout(
            [[length_hist, q_hist], [bc_counts]],
            sizing_mode="stretch_width"))

    section = report_doc.add_section()
    section.markdown('''
### Genome coverage
Plots below indicate depth of coverage from data used within the Artic analysis
coloured by amplicon pool. For adequate variant calling depth should be at least
30X in any region. Pool-1 reads are shown in light-blue, Pool-2 reads are dark grey
(a similarly for forward and reverse reads respectively).
''')

    # depth summary by amplicon pool
    df = read_files(args.depths)
    plots_pool = list()
    plots_orient = list()
    depth_lim = 100
    for sample in sorted(df['sample_name'].unique()):
        bc = df['sample_name'] == sample
        depth = df[bc].groupby('pos')['depth'].sum()
        depth_thresh = 100*(depth >= depth_lim).sum() / len(depth)

        # primer set plot
        pset = df['primer_set']
        xs = [df.loc[(pset == i) & bc]['pos'] for i in (1,2)]
        ys = [df.loc[(pset == i) & bc]['depth'] for i in (1,2)]

        print("------------", sample)
        print(df)
        plot = points.points(
            xs, ys, colors=[Colors.light_cornflower_blue, Colors.feldgrau],
            title="{}: {:.0f}X, {:.1f}% > {}X".format(
                sample, depth.mean(), depth_thresh, depth_lim),
            height=200, width=400,
            x_axis_label='position', y_axis_label='depth',
            ylim=(0,300))
        plots_pool.append(plot)

        # fwd/rev
        data = df[bc].groupby('pos').sum().reset_index()
        xs = [data['pos'], data['pos']]
        ys = [data['depth_fwd'], data['depth_rev']]
        plot = points.points(
            xs, ys, colors=[Colors.light_cornflower_blue, Colors.feldgrau],
            title="{}: {:.0f}X, {:.1f}% > {}X".format(
                sample, depth.mean(), depth_thresh, depth_lim),
            height=200, width=400,
            x_axis_label='position', y_axis_label='depth',
            ylim=(0,300))
        plots_orient.append(plot)

    tab1 = Panel(
        child=gridplot(plots_pool, ncols=3), title="By amplicon pool")
    tab2 = Panel(
        child=gridplot(plots_orient, ncols=3), title="By read orientation")
    cover_panel = Tabs(tabs=[tab1, tab2])
    section.plot(cover_panel)

    # canned VCF stats report component
    section = report_doc.add_section()
    bcfstats.full_report(args.bcftools_stats, report=section)

    section = report_doc.add_section()
    section.markdown('''
### NextClade analysis
The following view is produced by the [nextclade](https://clades.nextstrain.org/) software.
''')
    with open(args.nextclade, encoding='utf8') as fh:
        nc = fh.read()
    nextclade = report.NextClade(nc)
    section.plot(nextclade)

    # Footer section
    section = report_doc.add_section()
    section.markdown('''
### About

**Oxford Nanopore Technologies products are not intended for use for health assessment
or to diagnose, treat, mitigate, cure or prevent any disease or condition.**

This report was produced using the [epi2me-labs/wf-artic](https://github.com/epi2me-labs/wf-artic).
The workflow can be run using `nextflow epi2me-labs/wf-artic --help`

---
''')

    # write report
    report_doc.write("summary_report.html")

if __name__ == "__main__":
    main()
