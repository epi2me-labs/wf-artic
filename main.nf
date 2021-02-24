#!/usr/bin/env extflow

nextflow.enable.dsl = 2

valid_schemes = ["V1", "V2", "V3", "V1200"]

def helpMessage(){
    log.info """
SARS-Cov-2 Artic Analysis Workflow

Usage:
    nextflow run main.nf [options]

Options:
    --fastq             DIR     Path to FASTQ directory (required)
    --samples           FILE    CSV file with columns named `barcode` and `sample_name`
                                (or simply a sample name for non-multiplexed data).
    --out_dir           DIR     Path for output (default: $params.out_dir)
    --medaka_model      STR     Medaka model name (default: $params.medaka_model)
    --min_len           INT     Minimum read length (default: set by scheme)
    --max_len           INT     Maximum read length (default: set by scheme)
    --scheme_version    STR     Primer scheme ($valid_schemes)
                                indicating correspondence between
                                barcodes and sample names. (default: $params.scheme_version)

Notes:
    If directories named "barcode*" are found under the `--fastq` directory the
    data is assumed to be multiplex and each barcode directory will be processed
    independently. If `.fastq(.gz)` files are found under the `--fastq` directory
    the sample is assumed to not be multiplexed. In this second case `--samples`
    should be a simple name rather than a CSV file.

    Minimum and maximum rad length filters are applied based on the amplicon scheme.
    These can be overridden using the `--min_len` and `--max_len` options.
"""
}

process preArticQC {
    label "artic"
    cpus 1
    input:
        tuple file(directory), val(sample_name) 
    output:
        file "seq_summary.txt"
    """
#!/usr/bin/env python
import glob
import itertools
import os
import pysam
import numpy as np

fastqs = glob.glob(os.path.join("$directory", "*.fastq*"))
reads = itertools.chain.from_iterable(
    pysam.FastxFile(fname) for fname in fastqs)

def mean_qual(quals):
    qual = np.fromiter(
        (ord(x) - 33 for x in quals),
        dtype=int, count=len(quals))
    mean_p = np.mean(np.power(10, qual / -10))
    return -10 * np.log10(mean_p)

with open("seq_summary.txt", "w") as fh:
    # names as in Guppy
    fh.write("sample_name\\tread_id\\tsequence_length_template\\tmean_qscore_template\\n")
    for read in reads:
        fh.write("\\t".join(str(x) for x in (
            "$sample_name", read.name, len(read.sequence), mean_qual(read.quality))))
        fh.write("\\n") 
    """
}


process runArtic {

    label "artic"
    cpus 2
    input:
        tuple file(directory), val(sample_name)
        file "primer_schemes"
    output:
        file("${sample_name}.consensus.fasta")
        tuple file("${sample_name}.pass.named.vcf.gz"), file("${sample_name}.pass.named.vcf.gz.tbi")
        file("${sample_name}.depth.txt")

    """
    echo "=========="
    echo $directory, $sample_name
    echo "=========="

    function mock_artic {
        # write an empty VCF
        echo "Mocking artic results"
    cat << EOF |
    ##fileformat=VCFv4.2
    ##source=Longshot v0.4.0
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
    EOF
        bgzip > $sample_name".pass.vcf.gz"
        touch $sample_name".consensus.fasta"
    }

    # name everything by the sample rather than barcode
    if [[ "$directory" != "$sample_name" ]]; then
        echo "Moving input: '$directory' to '$sample_name'"
        mv $directory $sample_name
    fi

    artic guppyplex --skip-quality-check \
        --min-length $params._min_len --max-length $params._max_len \
        --directory $sample_name --prefix $sample_name \
        && echo " - artic guppyplex finished"
    # the output of the above will be...
    READFILE=$sample_name"_"$sample_name".fastq"

    artic minion --medaka --normalise 200 --threads $task.cpus \
        --read-file \$READFILE \
        --medaka-model $params.medaka_model \
        --scheme-directory primer_schemes \
        $params.full_scheme_name $sample_name \
        || mock_artic

    zcat $sample_name".pass.vcf.gz" | sed "s/SAMPLE/$sample_name/" | bgzip > $sample_name".pass.named.vcf.gz"
    bcftools index -t $sample_name".pass.named.vcf.gz"

    # rename the consensus sequence
    sed -i "s/^>.*/>$sample_name/" $sample_name".consensus.fasta"

    # calculate depth stats. Final output is single file annotated with primer set and sample name
    for i in 1 2; do
        bam=$sample_name".primertrimmed.\$i.sorted.bam"
        samtools view -r \$i -b $sample_name".primertrimmed.rg.sorted.bam" > \$bam
        samtools index \$bam
        stats_from_bam \$bam > \$bam".stats" || echo "stats_from_bam failed, probably no alignments"
        coverage_from_bam -s 20 -p \$bam \$bam
        # TODO: we're assuming a single reference sequence here
        awk 'BEGIN{OFS="\t"}{if(NR==1){print \$0, "sample_name", "primer_set"}else{print \$0, "$sample_name", '\$i'}}' *\${bam}*".depth.txt" > $sample_name".depth.\$i.txt"
        rm -rf \$bam \$bam.bai
    done
    cat $sample_name".depth.1.txt" <(tail -n+2 $sample_name".depth.2.txt") > $sample_name".depth.txt"
    """
}


process report {
    label "artic"
    cpus 1
    input:
        file "depths_*.txt"
        file "read_summary_*.txt"
        file "nextclade.json"
    output:
        file "summary_report.html"
    """
#!/usr/bin/env python

import glob
import numpy as np
import pandas as pd

from bokeh.layouts import gridplot, layout
from bokeh.models import Panel, Tabs
import aplanat
from aplanat import annot, bars, gridplot, hist, lines, points, report

report_doc = report.HTMLReport(
    "SARS-CoV-2 ARTIC Sequencing report",
    "Results generated through the wf-artic nextflow workflow provided by Oxford Nanopore Technologies")

def read_files(pattern):
    dfs = list()
    for fname in sorted(glob.glob(pattern)):
        dfs.append(pd.read_csv(fname, sep="\\t"))
    return pd.concat(dfs)

report_doc.markdown('''
### Read Quality control
This section displays basic QC metrics indicating read data quality.
''')

min_read_length = $params._min_len
max_read_length = $params._max_len
np_blue = '#0084A9'
np_dark_grey = '#455560'
np_light_blue = '#90C6E7'


# read length summary
seq_summary = read_files("read_summary_*.txt")
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
    length_hist, min_read_length,
    label="Min: {}".format(min_read_length), text_baseline='bottom', color='grey')
length_hist = annot.marker_vline(
    length_hist, max_read_length,
    label="Max: {}".format(max_read_length), text_baseline='top')
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
    (seq_summary['sequence_length_template'] > min_read_length)
    & (seq_summary['sequence_length_template'] < max_read_length)]
barcode_counts = (
    pd.DataFrame(good_reads['sample_name'].value_counts())
    .sort_index()
    .reset_index()
    .rename(
        columns={'index':'sample', 'sample_name':'count'})
    )

bc_counts = bars.simple_bar(
    barcode_counts['sample'].astype(str), barcode_counts['count'], colors=[np_blue]*len(barcode_counts),
    title='Number of reads per barcode (filtered by {} < length < {})'.format(min_read_length, max_read_length),
    plot_width=None
)
bc_counts.xaxis.major_label_orientation = 3.14/2
report_doc.plot(
    layout(
        [[length_hist, q_hist], [bc_counts]],
        sizing_mode="stretch_width"))

# depth summary by amplicon pool
df = read_files("depths_*.txt")
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

# nextclade data
with open("nextclade.json", encoding='utf8') as fh:
    nc = fh.read()
nextclade = report.NextClade(nc)
report_doc.markdown('''
### NextClade analysis
The following view is produced by the [nextclade](https://clades.nextstrain.org/) software.
''')
report_doc.plot(nextclade)

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
    """
}


process allConsensus {
    label "artic"
    cpus 1
    input:
        file "*"
    output:
        file "all_consensus.fasta"
    """
    ls *.consensus.fasta | xargs cat > all_consensus.fasta
    """
}


process allVariants {
    label "artic"
    cpus 1
    input:
        tuple file(vcfs), file(tbis)
    output:
        file "all_variants.vcf.gz"
        file "all_variants.vcf.gz.tbi"
    """
    bcftools merge -o all_variants.vcf.gz -O z *.vcf.gz
    bcftools index -t all_variants.vcf.gz 
    """
}


process nextclade {
    label "artic"
    cpus 1
    input:
        file "consensus.fasta"
        file "reference.fasta"
        file scheme_bed
    output:
        file "nextclade.json"
    """
    scheme_to_nextclade.py $scheme_bed reference.fasta primers.csv
    nextclade \
        --input-fasta consensus.fasta --input-pcr-primers primers.csv \
        --output-json nextclade.json
    """
}


// See https://github.com/nextflow-io/nextflow/issues/1636
// This is the only way to publish files from a workflow whilst
// decoupling the publish from the process steps.
process output {
    // publish inputs to output directory

    publishDir "${params.out_dir}", mode: 'copy', pattern: "*"
    input:
        file fname
    output:
        file fname
    """
    echo "Writing output files"
    """
}


// workflow module
workflow pipeline {
    take:
        samples
        scheme_directory
        reference
        primers
    main:
        read_summaries = preArticQC(samples)
        runArtic(samples, scheme_directory)
        // collate consensus and variants
        all_consensus = allConsensus(runArtic.out[0].collect())
        tmp = runArtic.out[1].toList().transpose().toList() // surely theres another way?
        all_variants = allVariants(tmp)
        // nextclade
        clades = nextclade(all_consensus, reference, primers)
        // report
        html_doc = report(
            runArtic.out[2].collect(), read_summaries.collect(), clades.collect())

        results = all_consensus.concat(all_variants, html_doc)
    emit:
        results
}

// entrypoint workflow
workflow {

    if (params.help) {
        helpMessage()
        exit 1
    }

    if (!params.fastq) {
        helpMessage()
        println("")
        println("`--fastq` is required")
        exit 1
    }

    if (!valid_schemes.any { it == params.scheme_version}) {
        println("`--scheme_version should be one of: $valid_schemes")
        exit 1
    }

    if (!params.min_len) {
        params.remove('min_len')
        if (params.scheme_version == "V1200") {
            params._min_len = 150
        } else {
            params._min_len = 400
        }
    } else {
        params._min_len = params.min_len
        params.remove('min_len')
    }
    if (!params.max_len) { 
        params.remove('max_len')
        if (params.scheme_version == "V1200") {
            params._max_len = 1200
        } else {
            params._max_len = 700
        }
    } else {
        params._min_len = params.min_len
        params.remove('min_len')
    }
    println("")
    println("Parameter summary")
    println("=================")
    params.each { it -> println("    $it.key: $it.value") }
    println("")

    params.full_scheme_name = params.scheme_name + "/" + params.scheme_version
    schemes = projectDir + '/data/primer_schemes' 
    scheme_directory = file(schemes, type: 'dir', checkIfExists:true)
    reference = file(
        "${scheme_directory}/${params.full_scheme_name}/${params.scheme_name}.reference.fasta",
        type:'file', checkIfExists:true)
    primers = file(
        "${scheme_directory}/${params.full_scheme_name}/${params.scheme_name}.scheme.bed",
        type:'file', checkIfExists:true)

    // resolve whether we have demultiplexed data or single sample
    barcode_dirs = file("$params.fastq/barcode*", type: 'dir', maxdepth: 1)
    not_barcoded = file("$params.fastq/*.fastq*", type: 'file', maxdepth: 1)
    if (barcode_dirs) {
        println("Found barcode directories")
        if (params.samples) {
            sample_sheet = Channel
                .fromPath(params.samples, checkIfExists: true)
                .splitCsv(header: true)
                .map { row -> tuple(row.barcode, row.sample_name) }
        } else {
            // just map directory name to self
            sample_sheet = Channel
                .fromPath(barcode_dirs)
                .filter(~/.*barcode[0-9]{1,3}$/)  // up to 192
                .map { path -> tuple(path.baseName, path.baseName) }
        }
        Channel
            .fromPath(barcode_dirs)
            .filter(~/.*barcode[0-9]{1,3}$/)  // up to 192
            .map { path -> tuple(path.baseName, path) }
            .join(sample_sheet)
            .map { barcode, path, sample -> tuple(path, sample) }
            .set{ samples }
    } else if (not_barcoded) {
        println("Found fastq files, assuming single sample")
        sample_sheet = Channel.of([$params.fastq, $params.samples])
        Channel
            .fromPath($params.fastq, type: 'dir', maxdepth:1)
            .join(sample_sheet)
            .set{ samples }
    }
    results = pipeline(samples, scheme_directory, reference, primers)
    output(results)
}
