#!/usr/bin/env python

import argparse
import glob
import os

import aplanat
from aplanat import annot, bars, gridplot, hist, lines, points, report
from bokeh.layouts import gridplot, layout
from bokeh.models import Panel, Tabs
import numpy as np
import pandas as pd

from parse import parse_bcftools_stats_multi

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

    report_doc.markdown('''
### Read Quality control
This section displays basic QC metrics indicating read data quality.
''')

    np_blue = '#0084A9'
    np_dark_grey = '#455560'
    np_light_blue = '#90C6E7'

    # read length summary
    seq_summary = read_files(args.summaries)
    total_bases = seq_summary['sequence_length_template'].sum()
    mean_length = total_bases / len(seq_summary)
    median_length = np.median(seq_summary['sequence_length_template'])
    datas = [seq_summary['sequence_length_template']]
    length_hist = hist.histogram(
        datas, colors=[np_blue], bins=100,
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
        datas, colors=[np_blue], bins=100,
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
        barcode_counts['sample'].astype(str), barcode_counts['count'], colors=[np_blue]*len(barcode_counts),
        title='Number of reads per barcode (filtered by {} < length < {})'.format(args.min_len, args.max_len),
        plot_width=None
    )
    bc_counts.xaxis.major_label_orientation = 3.14/2
    report_doc.plot(
        layout(
            [[length_hist, q_hist], [bc_counts]],
            sizing_mode="stretch_width"))

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

        plot = points.points(
            xs, ys, colors=[np_light_blue, np_dark_grey],
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
            xs, ys, colors=[np_light_blue, np_dark_grey],
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

    report_doc.markdown('''
### Genome coverage
Plots below indicate depth of coverage from data used within the Artic analysis
coloured by amplicon pool. For adequate variant calling depth should be at least
30X in any region. Pool-1 reads are shown in light-blue, Pool-2 reads are dark grey
(a similarly for forward and reverse reads respectively).
''')
    report_doc.plot(cover_panel)

    report_doc.markdown('''
### Variant call summaries

The following tables and figures are derived from the output of `bcftools stats`.
''')
    vcf_tables = parse_bcftools_stats_multi(args.bcftools_stats)
    report_doc.markdown("""
**Variant counts:**
""")
    df = vcf_tables['SN'] \
        .drop(columns='samples').set_index('sample').transpose()
    report_doc.table(df, index=True)
    report_doc.markdown("**Transitions and tranversions:**")
    df = vcf_tables['TSTV'] \
        .set_index('sample').transpose()
    report_doc.table(df, index=True)
    report_doc.markdown("""
**Substitution types**

Base substitutions aggregated across all samples (symmetrised by pairing)
""")

    sim_sub = {
        'G>A': 'C>T', 'G>C': 'C>G', 'G>T': 'C>A',
        'T>A': 'A>T', 'T>C': 'A>G', 'T>G': 'A>C'}
    def canon_sub(sub):
        b1 = sub[0]
        if b1 not in {'A', 'C'}:
            return canon_sub(sim_sub[sub])
        else:
            return b1, sub[2]

    df = vcf_tables['ST']
    df['canon_sub'] = df['type'].apply(canon_sub)
    df['original'] = df['canon_sub'].apply(lambda x: x[0])
    df['substitution'] = df['canon_sub'].apply(lambda x: x[1])
    df['count'] = df['count'].astype(int)
    df = df[['original', 'substitution', 'count']] \
        .groupby(['original', 'substitution']) \
        .agg(count=pd.NamedAgg(column='count', aggfunc='sum')) \
        .reset_index()

    from bokeh.models import ColorBar, LinearColorMapper
    from bokeh.palettes import Blues9
    from bokeh.plotting import figure
    colors = Blues9[::-1]
    mapper = LinearColorMapper(
        palette=colors, low=min(df['count']), high=max(df['count']))
    p = figure(
        y_range=['C', 'A'], x_range=['A', 'C', 'G', 'T'],
        x_axis_location="above",
        x_axis_label='alternative base',
        y_axis_label='reference base',
        tools="save", toolbar_location='below',
        output_backend="webgl",
        height=225, width=300,
        tooltips=[('sub', '@original>@substitution'), ('count', '@count')])
    p.grid.grid_line_color = None
    p.axis.axis_line_color = None
    p.axis.major_tick_line_color = None
    p.rect(
        source=df, y="original", x="substitution", width=1, height=1,
        fill_color={'field': 'count', 'transform': mapper},
        line_color=None)
    color_bar = ColorBar(
        title='', color_mapper=mapper, label_standoff=10,
        location=(0, 0))
    #p.add_layout(color_bar, 'right')

    report_doc.plot(p)

    report_doc.markdown("""
**Indel lengths**

Insertion and deletion lengths aggregated across all samples.
""")
    try:
        df = vcf_tables['IDD']
    except KeyError as e:
        # If there are no indels, bcftools doesn't contain the table
        report_doc.markdown("*No indels to report.*")
    else:
        df['nlength'] = df['length (deletions negative)'].astype(int)
        df['count'] = df['number of sites'].astype(int)
        # pad just to pull out axes by a minimum
        pad = pd.DataFrame({'nlength':[-10,+10], 'count':[0,0]})
        counts = df.groupby('nlength') \
            .agg(count=pd.NamedAgg(column='count', aggfunc='sum')) \
            .reset_index().append(pad)
        plot = hist.histogram(
            [counts['nlength']], weights=[counts['count']],
            colors = [np_light_blue], binwidth=1,
            title='Insertion and deletion variant lengths',
            x_axis_label='Length / bases (deletions negative)',
            y_axis_label='Count')
        #plot.xaxis.formatter = NumeralTickFormatter(format="0,0")
        report_doc.plot(plot)

    report_doc.markdown('''
### NextClade analysis
The following view is produced by the [nextclade](https://clades.nextstrain.org/) software.
''')
    with open(args.nextclade, encoding='utf8') as fh:
        nc = fh.read()
    nextclade = report.NextClade(nc)
    report_doc.plot(nextclade)

    # Footer section
    report_doc.markdown('''
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
