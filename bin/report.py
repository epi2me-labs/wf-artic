#!/usr/bin/env python
"""Create report file."""

import argparse

from aplanat import annot, bars, hist, lines, report
from aplanat.components import bcfstats, nextclade
from aplanat.util import Colors
from bokeh.layouts import gridplot, layout
from bokeh.models import Panel, Range1d, Tabs
import conda_versions
import numpy as np
import pandas as pd


def load_params(path):
    """Load the parameters csv into a dataframe."""
    params_cols = ['Name', 'Value']
    df_params = pd.read_csv(path, names=params_cols)

    for _, row in df_params.iterrows():
        name = row['Name']
        value = row['Value']

        if not name.startswith('_'):
            continue

        initial = name.lstrip('_')
        matches = df_params.loc[df_params['Name'] == initial]
        if matches.empty:
            df_new = pd.DataFrame([[initial, value]], columns=params_cols)
            df_params = df_params.append(df_new, ignore_index=True)
        else:
            matches['Value'] = value

        df_params = df_params[df_params['Name'] != name]

    df_params.fillna('N/A', inplace=True)
    return df_params


def read_files(summaries, sep='\t'):
    """Read a set of files and join to single dataframe."""
    dfs = list()
    for fname in sorted(summaries):
        dfs.append(pd.read_csv(fname, sep=sep))
    return pd.concat(dfs)


def main():
    """Run entry point."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "status",
        help="artic status file")
    parser.add_argument(
        "output",
        help="Report output filename")
    parser.add_argument(
        "--nextclade",
        help="nextclade json output file")
    parser.add_argument(
        "--pangolin",
        help="pangolin CSV output file")
    parser.add_argument(
        "--depths", nargs='+', required=True,
        help="Depth summary files")
    parser.add_argument(
        "--summaries", nargs='+', required=True,
        help="Sequencing summary files")
    parser.add_argument(
        "--bcftools_stats", nargs='+', required=True,
        help="Outputs from bcftools stats")
    parser.add_argument(
        "--genotypes", nargs='+', required=False,
        help="Genotyping summary files")
    parser.add_argument(
        "--min_cover", default=20, type=int,
        help="Minimum locus coverage for variant call.")
    parser.add_argument(
        "--min_len", default=300, type=int,
        help="Minimum read length")
    parser.add_argument(
        "--max_len", default=700, type=int,
        help="Maximum read length")
    parser.add_argument(
        "--report_depth", default=100, type=int,
        help=(
            "Depth at which to provide a coverage statistics, "
            "e.g. 76% of genome covered at `report_depth`"))
    parser.add_argument(
        "--hide_coverage", action="store_true",
        help="Do not display coverage plots in report.")
    parser.add_argument(
        "--hide_variants", action="store_true",
        help="Do not display variant summary in report.")
    parser.add_argument(
        "--revision", default='unknown',
        help="git branch/tag of the executed workflow")
    parser.add_argument(
        "--commit", default='unknown',
        help="git commit of the executed workflow")
    parser.add_argument(
        "--params", default=None,
        help="A csv containing the parameter key/values")
    args = parser.parse_args()

    report_doc = report.WFReport(
        "SARS-CoV-2 ARTIC Sequencing report", "wf-artic",
        revision=args.revision, commit=args.commit)

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
        datas, colors=[Colors.cerulean], binwidth=50,
        title="Read length distribution.",
        x_axis_label='Read Length / bases',
        y_axis_label='Number of reads',
        xlim=(0, 2000))
    length_hist = annot.marker_vline(
        length_hist, args.min_len,
        label="Min: {}".format(args.min_len), text_baseline='bottom',
        color='grey')
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
            columns={'index': 'sample', 'sample_name': 'count'})
    )

    bc_counts = bars.simple_bar(
        barcode_counts['sample'].astype(str), barcode_counts['count'],
        colors=[Colors.cerulean]*len(barcode_counts),
        title=(
            'Number of reads per barcode '
            '(filtered by {} < length < {})'.format(
                args.min_len, args.max_len)),
        plot_width=None
    )
    bc_counts.xaxis.major_label_orientation = 3.14/2
    section.plot(
        layout(
            [[length_hist, q_hist], [bc_counts]],
            sizing_mode="stretch_width"))

    section = report_doc.add_section()
    section.markdown("""
### Artic Analysis status

The panel below lists samples which failed to produce
results from the primary ARTIC analysis. Samples not listed here were analysed
successfully, but may still contain inconclusive or invalid results. See the
following sections for further indications of failed or inconclusive results.
""")
    status = pd.read_csv(args.status, sep='\t')
    failed = status.loc[status['pass'] == 0]
    if len(failed) == 0:
        fail_list = "All samples analysed successfully"
    else:
        fail_list = failed['sample'].str.cat(sep=', ')
    section.markdown("""
```{}```
""".format(fail_list))
    fail_percentage = int(100 * len(failed) / len(status))
    classes = ['Success', 'Analysis Failed']
    values = [100 - fail_percentage, fail_percentage]
    colors = ['#54B8B1', '#EF4135']
    plot = bars.single_hbar(
        values, classes, colors,
        title="Completed analyses",
        x_axis_label="%age Samples")
    plot.x_range = Range1d(0, 140)
    section.plot(plot)

    if not args.hide_coverage:
        section = report_doc.add_section()
        section.markdown('''
### Genome coverage

Plots below indicate depth of coverage from data used
within the Artic analysis coloured by amplicon pool.  Variant filtering during
the ARTIC analysis mandates a minimum coverage of at least {}X at
variant/genotyping loci for a call to be made.

***NB: To better display all possible data, the depth axes of the plots below
are not tied between plots for different samples. Care should be taken in
comparing depth across samples.***
'''.format(args.min_cover))

        # depth summary by amplicon pool
        df = read_files(args.depths)
        plots_pool = list()
        plots_orient = list()
        plots_combined = list()
        depth_lim = args.report_depth
        for sample in sorted(df['sample_name'].unique()):
            bc = df['sample_name'] == sample
            depth = df[bc].groupby('pos').sum().reset_index()
            depth_thresh = \
                100*(depth['depth'] >= depth_lim).sum() / len(depth['depth'])
            depth_mean = depth['depth'].mean()

            # total depth plot
            # plot line just to get aplanat niceities
            p = lines.line(
                [depth['pos']], [depth['depth']], colors=[Colors.cerulean],
                title="{}: {:.0f}X, {:.1f}% > {}X".format(
                    sample, depth_mean, depth_thresh, depth_lim),
                height=250, width=400,
                x_axis_label='position', y_axis_label='depth',
            )
            p.varea(
                x=depth['pos'], y1=0.1, y2=depth['depth'],
                fill_color=Colors.cerulean)
            plots_combined.append(p)

            # fwd/rev
            xs = [depth['pos'], depth['pos']]
            ys = [depth['depth_fwd'], depth['depth_rev']]
            names = ['fwd', 'rev']
            colors = [Colors.dark_gray, Colors.verdigris]

            p = lines.line(
                xs, ys, colors=colors, names=names,
                title="{}: {:.0f}X, {:.1f}% > {}X".format(
                    sample, depth_mean, depth_thresh, depth_lim),
                height=250, width=400,
                x_axis_label='position', y_axis_label='depth')
            for x, y, name, color in zip(xs, ys, names, colors):
                p.varea(
                    x=x, y1=0, y2=y, legend_label=name,
                    fill_color=color, alpha=0.7,
                    muted_color=color, muted_alpha=0.2)
            p.legend.click_policy = 'mute'
            plots_orient.append(p)

            # primer set plot
            pset = df['primer_set']
            xs = [df.loc[(pset == i) & bc]['pos'] for i in (1, 2)]
            ys = [df.loc[(pset == i) & bc]['depth'] for i in (1, 2)]
            names = ['pool-1', 'pool-2']
            colors = [Colors.light_cornflower_blue, Colors.feldgrau]

            p = lines.line(
                xs, ys, colors=colors, names=names,
                title="{}: {:.0f}X, {:.1f}% > {}X".format(
                    sample, depth_mean, depth_thresh, depth_lim),
                height=250, width=400,
                x_axis_label='position', y_axis_label='depth')
            for x, y, name, color in zip(xs, ys, names, colors):
                p.varea(
                    x=x, y1=0, y2=y, legend_label=name,
                    fill_color=color, alpha=0.7,
                    muted_color=color, muted_alpha=0.2)
            p.legend.click_policy = 'mute'
            plots_pool.append(p)

        tab1 = Panel(
            child=gridplot(plots_combined, ncols=3), title="Coverage Plot")
        tab2 = Panel(
            child=gridplot(plots_pool, ncols=3), title="By amplicon pool")
        tab3 = Panel(
            child=gridplot(plots_orient, ncols=3), title="By read orientation")
        cover_panel = Tabs(tabs=[tab1, tab2, tab3])
        section.plot(cover_panel)

    # canned VCF stats report component
    if not args.hide_variants:
        section = report_doc.add_section()
        bcfstats.full_report(args.bcftools_stats, report=section)

    # NextClade analysis
    if args.nextclade is not None:
        section = report_doc.add_section(
            section=nextclade.NextClade(args.nextclade))
        section.markdown(
            "*Note: For targeted sequencing, such as SpikeSeq, Nextclade "
            "may report 'Missing data' QC fails. This is expected and not "
            "a concern provided the regions of interest are not reported "
            "as missing.*")

    # Pangolin analysis
    if args.pangolin is not None:
        section = report_doc.add_section()
        section.markdown('''
### Lineage

The table below reports the lineage of each sample as calculated by
[pangolin](https://github.com/cov-lineages/pangolin).

''')
        section.table(pd.read_csv(args.pangolin), index=False)

    # Genotyping
    if args.genotypes is not None:
        section = report_doc.add_section()
        section.markdown('''
### Genotyping

The table below lists whether candidate variants were determined to exist
within each sample.

The ARTIC workflow pre-filters (removes) candidate variants according to the
criteria `variant_score < 20` and `coverage < 20`. The table draws attention to
reference calls of low coverage (<20 reads) which may therefore be inaccurate.
''')
        df = read_files(args.genotypes, sep=',')
        df = df[[
            'Sample', 'Result', 'Date Tested', 'Lab ID', 'testKit',
            'CH1-Target', 'CH1-Result', 'CH1-Conf']]
        df = df.sort_values(by=['Sample', 'CH1-Target'], ascending=True)
        section.table(df, index=False)

    section = report_doc.add_section()
    section.markdown('''
### Software versions

The table below highlights versions of key software used within the analysis.
''')
    req = [
        'np-artic', 'medaka', 'minimap2', 'bcftools', 'nextclade-cli',
        'samtools', 'bcftools']
    versions = conda_versions.scrape_data(
        as_dataframe=True, include=req)
    section.table(versions[['Name', 'Version', 'Build']], index=False)

    # Params reporting
    section = report_doc.add_section()
    section.markdown('''
### Workflow parameters

The table below highlights values of the main parameters used in this analysis.
''')
    df_params = load_params(args.params)
    section.table(df_params, index=False)

    # write report
    report_doc.write(args.output)


if __name__ == "__main__":
    main()
