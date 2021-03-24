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
    --report_depth      INT     Min. depth for percentage coverage (default: $params.report_depth)
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


process checkSampleSheet {
    label "artic"
    cpus 1
    input:
        file "sample_sheet.txt"
    output:
        file "samples.txt"
    """
#!/usr/bin/env python
import pandas as pd
try:
    samples = pd.read_csv("sample_sheet.txt", sep=None)
    if 'barcode' not in samples.columns or 'sample_name' not in samples.columns:
        raise IOError()
except Exception:
    raise IOError(
        "Could not parse sample sheet, it must contain two columns "
        "named 'barcode' and 'sample_name'.")
# check duplicates
dup_bc = samples['barcode'].duplicated()
dup_sample = samples['sample_name'].duplicated()
if any(dup_bc) or any(dup_sample):
    raise IOError(
        "Sample sheet contains duplicate values.")
samples.to_csv("samples.txt", sep=",", index=False)
    """
}

process preArticQC {
    label "artic"
    cpus 1
    input:
        tuple file(directory), val(sample_name) 
    output:
        file "${sample_name}.stats"
    """
    seq_summary.py $directory ${sample_name}.stats --sample_name $sample_name
    """
}


process runArtic {
    label "artic"
    cpus 2
    input:
        tuple file(directory), val(sample_name)
        file "primer_schemes"
    output:
        file "${sample_name}.consensus.fasta"
        tuple file("${sample_name}.pass.named.vcf.gz"), file("${sample_name}.pass.named.vcf.gz.tbi")
        file "${sample_name}.depth.txt"
        file "${sample_name}.pass.named.stats"
    """
    run_artic.sh \
        ${sample_name} ${directory} ${params._min_len} ${params._max_len} \
        ${params.medaka_model} ${params.full_scheme_name} \
        ${task.cpus}
    bcftools stats ${sample_name}.pass.named.vcf.gz > ${sample_name}.pass.named.stats 
    """
}


process report {
    label "artic"
    cpus 1
    input:
        file "depth_stats/*"
        file "read_stats/*"
        file "nextclade.json"
        file "vcf_stats/*"
        file "consensus_status.txt"
    output:
        file "wf-artic-report.html"
    """
    report.py nextclade.json consensus_status.txt wf-artic-report.html \
        --min_len $params._min_len --max_len $params._max_len --report_depth \
        $params.report_depth --depths depth_stats/* --summaries read_stats/* \
        --bcftools_stats vcf_stats/*
    """
}


process allConsensus {
    label "artic"
    cpus 1
    input:
        file "*"
    output:
        file "all_consensus.fasta"
        file "consensus_status.txt"
    """
    ls *.consensus.fasta | xargs cat > all_consensus.fasta
    grep "^>" all_consensus.fasta \
        | awk 'BEGIN{OFS="\\t"; print "sample\\tpass"}{print substr(\$1, 2), \$2!="Artic-Fail"}' \
        >> consensus_status.txt
    """
}


process allVariants {
    label "artic"
    cpus 1
    input:
        tuple file(vcfs), file(tbis)
    output:
        tuple file("all_variants.vcf.gz"), file("all_variants.vcf.gz.tbi")
    """
    if [[ \$(ls *.vcf.gz | wc -l) == "1" ]]; then
        mv *.vcf.gz all_variants.vcf.gz
        mv *.vcf.gz.tbi all_variants.vcf.gz.tbi
    else 
        bcftools merge -o all_variants.vcf.gz -O z *.vcf.gz
        bcftools index -t all_variants.vcf.gz
    fi
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
        clades = nextclade(all_consensus[0], reference, primers)
        // report
        html_doc = report(
            runArtic.out[2].collect(),
            read_summaries.collect(), 
            clades.collect(), 
            runArtic.out[3].collect(),
            all_consensus[1])
        results = all_consensus[0].concat(all_consensus[1], all_variants[0].flatten(), html_doc)
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
    // check sample sheet
    sample_sheet = null
    if (params.samples) {
        sample_sheet = Channel.fromPath(params.samples, checkIfExists: true)
        sample_sheet = checkSampleSheet(sample_sheet)
            .splitCsv(header: true)
            .map { row -> tuple(row.barcode, row.sample_name) }
    }

    // resolve whether we have demultiplexed data or single sample
    not_barcoded = file("$params.fastq/*.fastq*", type: 'file', maxdepth: 1)
    // remove empty barcode_dirs
    barcode_dirs = file("$params.fastq/barcode*", type: 'dir', maxdepth: 1)
    valid_barcode_dirs = []
    invalid_barcode_dirs = []
    for (d in barcode_dirs) {
        if(!file("${d}/*.fastq*", type:'file', checkIfExists:true)) {
            invalid_barcode_dirs << d
        } else {
            valid_barcode_dirs << d
        }
    }
    if (barcode_dirs) {
        println("Found barcode directories")
        if (invalid_barcode_dirs.size() > 0) {
            println("Some barcode directories did not contain .fastq(.gz) files:")
            for (d in invalid_barcode_dirs) {
                println("- ${d}")
            }
        }
        if (!sample_sheet) {
            // just map directory name to self
            sample_sheet = Channel
                .fromPath(valid_barcode_dirs)
                .filter(~/.*barcode[0-9]{1,3}$/)  // up to 192
                .map { path -> tuple(path.baseName, path.baseName) }
        }
        Channel
            .fromPath(valid_barcode_dirs)
            .filter(~/.*barcode[0-9]{1,3}$/)  // up to 192
            .map { path -> tuple(path.baseName, path) }
            .join(sample_sheet)
            .map { barcode, path, sample -> tuple(path, sample) }
            .set{ samples }
    } else if (not_barcoded) {
        println("Found fastq files, assuming single sample")
        sample = (params.samples == null) ? "unknown" : params.samples
        Channel
            .fromPath(params.fastq, type: 'dir', maxDepth:1)
            .map { path -> tuple(path, sample) }
            .set{ samples }
    }
    results = pipeline(samples, scheme_directory, reference, primers)
    output(results)
}
