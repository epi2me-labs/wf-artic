#!/usr/bin/env extflow

nextflow.enable.dsl = 2

params.help = ""
params.out_dir = "output"
params.prefix = "artic"
params.min_len = "300"
params.max_len = "700"
params.medaka_model = "r941_min_high_g360"
params.scheme_name = "SARS-CoV-2"
params.scheme_version = "V3"


if(params.help) {
    log.info ''
    log.info 'Workflow template'
    log.info ''
    log.info 'Usage: '
    log.info '    nextflow run workflow.nf [options]'
    log.info ''
    log.info 'Script Options: '
    log.info '    --fastq             DIR     Path to FASTQ directory'
    log.info '    --out_dir           DIR     Path for output'
    log.info '    --prefix            STR     Output filename prefix'
    log.info '    --medaka_model      STR     Medaka model name'
    log.info '    --min_length        INT     Minimum read length'
    log.info '    --max_length        INT     Maximum read length'
    log.info '    --samples           FILE    CSV file with columns named `barcode` and `sample_name`'
    log.info '    --scheme_version    STR     Primer scheme (V1, V2, V3, V1200)'
    log.info '                                indicating correspondence between'
    log.info '                                barcodes and sample names.'
    log.info ''

    return
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
    ln -s $directory $sample_name
    artic guppyplex --skip-quality-check \
        --min-length $params.min_len --max-length $params.max_len \
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
        coverage_from_bam -s 10 -p \$bam \$bam
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

from bokeh.layouts import gridplot
import aplanat
from aplanat import annot, bars, gridplot, hist, lines, report

report_doc = report.HTMLReport(
    "SARS-CoV-2 ARTIC Sequencing report",
    "Results generated through the wf-artic nextflow workflow provided by Oxford Nanopore Technologies")

def read_files(pattern):
    dfs = list()
    for fname in glob.glob(pattern):
        dfs.append(pd.read_csv(fname, sep="\\t"))
    return pd.concat(dfs)

report_doc.markdown('''
#### Read Quality control
This section displays basic QC metrics indicating read data quality.
''')

min_read_length = $params.min_len
max_read_length = $params.max_len

# read length summary
seq_summary = read_files("read_summary_*.txt")
total_bases = seq_summary['sequence_length_template'].sum()
mean_length = total_bases / len(seq_summary)
median_length = np.median(seq_summary['sequence_length_template'])
datas = [seq_summary['sequence_length_template']]
length_hist = hist.histogram(
    datas, bins=400,
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
    datas, bins=400,
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
    barcode_counts['sample'].astype(str), barcode_counts['count'],
    title='Number of reads per barcode (filtered by {} < length < {})'.format(min_read_length, max_read_length))
bc_counts.xaxis.major_label_orientation = 3.14/2
report_doc.plot(gridplot([length_hist, q_hist, bc_counts], ncols=2))

# depth summary by amplicon pool
df = read_files("depths_*.txt")
plots = list()
depth_lim = 100
for sample in df['sample_name'].unique():
    pset = df['primer_set']
    bc = df['sample_name'] == sample
    xs = [df.loc[(pset == i) & bc]['pos'] for i in (1,2)]
    ys = [df.loc[(pset == i) & bc]['depth'] for i in (1,2)]
    
    depth = df[bc].groupby('pos')['depth'].sum()
    depth_thresh = 100*(depth >= depth_lim).sum() / len(depth)

    plot = lines.line(
        xs, ys, colors=['blue', 'red'],
        title="{}: {:.0f}X, {:.1f}% > {}X".format(
            sample, depth.mean(), depth_thresh, depth_lim),
        height=200, width=400,
        x_axis_label='position', y_axis_label='depth',
        ylim=(0,300))
    plots.append(plot)
report_doc.markdown('''
#### Genome coverage
Plots below indicate depth of coverage from data used within the Artic analysis
coloured by amplicon pool. For adequate variant calling depth should be at least
30X in any region.
''')
report_doc.plot(gridplot(plots, ncols=3))

# nextclade data
with open("nextclade.json", encoding='utf8') as fh:
    nc = fh.read()
nextclade = report.NextClade(nc)
report_doc.markdown('''
### NextClade analysis
The following view is produced by the [nextclade](https://clades.nextstrain.org/) software.
''')
report_doc.plot(nextclade)

# write report
report_doc.write("summary_report.html")
    """
}


process allConsensus {
    label "artic"
    cpus 1
    input:
        file "sample_*.fastq"
    output:
        file "all_consensus.fasta"
    """
    cat sample_*.fastq > all_consensus.fasta
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
    label "nextclade"
    cpus 1
    input:
        file "consensus.fasta"
    output:
        file "nextclade.json"
    """
    nextclade --input-fasta 'consensus.fasta' --output-json 'nextclade.json'
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
    main:
        read_summaries = preArticQC(samples)
        runArtic(samples, scheme_directory)
        // collate consensus and variants
        all_consensus = allConsensus(runArtic.out[0].collect())
        tmp = runArtic.out[1].toList().transpose().toList() // surely theres another way?
        all_variants = allVariants(tmp)
        // nextclade
        clades = nextclade(all_consensus)
        // report
        html_doc = report(
            runArtic.out[2].collect(), read_summaries.collect(), clades.collect())
        

        results = all_consensus.concat(all_variants, html_doc)
    emit:
        results
}

// entrypoint workflow
workflow {

    params.full_scheme_name = params.scheme_name + "/" + params.scheme_version
    schemes = projectDir + '/data/primer_schemes' 
    scheme_directory = file(schemes, type: 'dir', checkIfExists:true)

    // resolve whether we have demultiplexed data or single sample
    barcode_dirs = file("$params.fastq/barcode*", type: 'dir', maxdepth: 1)
    not_barcoded = file("$params.fastq/*.fastq", type: 'file', maxdepth: 1)
    if (barcode_dirs) {
        println("Found barcode directories")
        sample_sheet = Channel
            .fromPath(params.samples, checkIfExists: true)
            .splitCsv(header: true)
            .map { row -> tuple(row.barcode, row.sample_name) }
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
    //samples.view()
    results = pipeline(samples, scheme_directory)
    output(results)
}
