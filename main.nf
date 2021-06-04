#!/usr/bin/env nextflow
import java.util.zip.GZIPInputStream;

nextflow.enable.dsl = 2

valid_schemes = ["SARS-CoV-2", "spike-seq"]
valid_scheme_versions = ["V1", "V2", "V3", "V1200"]

if (params.scheme_name == "spike-seq") {
    valid_scheme_versions = ["V1"]
}

def helpMessage(){
    log.info """
SARS-Cov-2 Artic Analysis Workflow

Usage:
    nextflow run main.nf [options]

Options:
    --fastq                     DIR     Path to FASTQ directory (required)
    --samples                   FILE    CSV file with columns named `barcode` and `sample_name`
                                        (or simply a sample name for non-multiplexed data).
    --out_dir                   DIR     Path for output (default: $params.out_dir)
    --medaka_model              STR     Medaka model name (default: $params.medaka_model)
    --min_len                   INT     Minimum read length (default: set by scheme)
    --max_len                   INT     Maximum read length (default: set by scheme)
    --scheme_name               STR     Scheme to use ($valid_schemes) (default: SARS-CoV-2)
    --scheme_version            STR     Primer scheme version ($valid_scheme_versions)
                                        (default: $params.scheme_version)
    --report_depth              INT     Min. depth for percentage coverage (default: $params.report_depth)
                                        (e.g. 89% genome covered at > `report_depth`)
                                        indicating correspondence between
    --genotype_variants         FILE    Report genotyping information for scheme's known variants of interest,
                                        optionally provide file path as argument.
    --detect_samples            BOOL    Automatically determine sample_id information from fastq
                                        header, replaces the --samples csv.
    --report_clade              BOOL    Show results of Nextclade analysis in report.
    --report_lineage            BOOL    Show results of Pangolin analysis in report.
    --report_coverage           BOOL    Show genome coverage traces in report.
    --report_variant_summary    BOOL    Show / hide variant information in report. (default: true)

Metadata:
    --timestamp                 STR     Timestamp for the genotyping report
    --lab_id                    STR     Lab_id for genotyping report
    --testkit                   STR     Testkit for genotyping report

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


def sampleFromHeader(fastq) {
    InputStream stream
    if (fastq.endsWith("gz")) {
        stream = new GZIPInputStream(
            new FileInputStream(fastq.toString()));
    } else {
        stream = new FileInputStream(fastq.toString());
    }

    Reader reader = new InputStreamReader(stream, "UTF-8");
    def line = reader.readLine()
    def found = line.findAll("sample_id=[A-Za-z0-9]+")
    reader.close()

    if (found) {
       return found[0].split('=')[1]
    }

    println("--detect_samples is on but no sample_id was found in fastq header")
    exit 1
}


process checkSampleSheet {
    label "artic"
    cpus 1
    input:
        file "sample_sheet.txt"
    output:
        file "samples.txt"
    """
    check_sample_sheet.py sample_sheet.txt samples.txt
    """
}


process copySchemeDir {
    label "artic"
    cpus 1
    input:
        path scheme_directory
    output:
        path "scheme_dir"
    """
    cp -RL $scheme_directory scheme_dir
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
        path "${sample_name}.consensus.fasta", emit: consensus
        tuple path("${sample_name}.pass.named.vcf.gz"), path("${sample_name}.pass.named.vcf.gz.tbi"), emit: pass_vcf
        tuple path("${sample_name}.merged.gvcf.named.vcf.gz"), path("${sample_name}.merged.gvcf.named.vcf.gz.tbi"), emit: merged_gvcf
        path "${sample_name}.depth.txt", emit: depth_stats
        path "${sample_name}.pass.named.stats", emit: vcf_stats
        tuple path("${sample_name}.sorted.bam"), path("${sample_name}.sorted.bam.bai"), emit: bam
    """
    run_artic.sh \
        ${sample_name} ${directory} ${params._min_len} ${params._max_len} \
        ${params.medaka_model} ${params.full_scheme_name} \
        ${task.cpus}
    bcftools stats ${sample_name}.pass.named.vcf.gz > ${sample_name}.pass.named.stats 
    """
}


process genotypeSummary {
    // Produce a genotype summary spreadsheet
    label "artic"
    cpus 1
    input:
        tuple file(vcf), file(tbi)
        tuple file(bam), file(bam_index)
        file "reference.vcf"
    output:
        file "*genotype.csv"
    script:
        def lab_id = params.lab_id ? "--lab_id ${params.lab_id}" : ""
        def testkit = params.testkit ? "--testkit ${params.testkit}" : ""
        def csvName = vcf.simpleName
    """
    genotype_summary.py -b $bam -v $vcf -d reference.vcf --sample $csvName $lab_id $testkit -o ${csvName}.genotype.csv
    """
}


process combineGenotypeSummaries {
    label "artic"
    cpus 1
    input:
        file "summary_*.csv"
    output:
        file "genotype_summary.csv"
    """
    combine_genotype_summaries.py -g *.csv -o genotype_summary.csv
    """
}


process report {
    label "artic"
    cpus 1
    input:
        file "depth_stats/*"
        file "read_stats/*"
        file "nextclade.json"
        file "pangolin.csv"
        file "genotypes/*"
        file "vcf_stats/*"
        file "consensus_status.txt"
    output:
        file "wf-artic-report.html"
    script:
    // when genotype_variants is false the channel contains a mock file
    def genotype = params.genotype_variants ? "--genotypes genotypes/*" : ""
    def nextclade = params.report_clade as Boolean ? "--nextclade nextclade.json" : ""
    def pangolin = params.report_lineage as Boolean ? "--pangolin pangolin.csv" : ""
    def coverage = params.report_coverage as Boolean ? "" : "--hide_coverage"
    def var_summary = params.report_variant_summary as Boolean ? "" : "--hide_variants"
    def paramsMap = params.toMapString()
        .replace("[", "")
        .replace("]", "")
        .replace(", ", "\n")
        .replace(":", ",");
    """
    echo "$pangolin"
    echo "$nextclade"
    echo "$paramsMap" > params.csv
    report.py \
        consensus_status.txt wf-artic-report.html \
        $pangolin $nextclade $coverage $var_summary \
        --revision $workflow.revision --params params.csv --commit $workflow.commitId \
        --min_len $params._min_len --max_len $params._max_len --report_depth \
        $params.report_depth --depths depth_stats/* --summaries read_stats/* \
        --bcftools_stats vcf_stats/* $genotype
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
        file reference
    output:
        tuple file("all_variants.vcf.gz"), file("all_variants.vcf.gz.tbi")
    """
    for vcf in \$(ls *.vcf.gz)
    do
        bcftools norm -c s -O z --fasta-ref $reference \$vcf > norm.\$vcf
        bcftools index -t norm.\$vcf
    done
    if [[ \$(ls norm.*.vcf.gz | wc -l) == "1" ]]; then
        mv norm.*.vcf.gz all_variants.vcf.gz
        mv norm.*.vcf.gz.tbi all_variants.vcf.gz.tbi
    else 
        bcftools merge -o all_variants.vcf.gz -O z norm.*.vcf.gz
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
    cp -L reference.fasta ref.fasta
    scheme_to_nextclade.py $scheme_bed ref.fasta primers.csv
    nextclade \
        --input-fasta consensus.fasta --input-pcr-primers primers.csv \
        --output-json nextclade.json --jobs 1
    """
}


process pangolin {
    label "pangolin"
    cpus 1
    input:
        file "consensus.fasta"
    output:
        file "lineage_report.csv"
    """
    pangolin consensus.fasta
    """
}


// See https://github.com/nextflow-io/nextflow/issues/1636
// This is the only way to publish files from a workflow whilst
// decoupling the publish from the process steps.
process output {
    // publish inputs to output directory
    label "artic"

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
        ref_variants
    main:
        combined_genotype_summary = null
        scheme_directory = copySchemeDir(scheme_directory)
        read_summaries = preArticQC(samples)
        runArtic(samples, scheme_directory)
        // collate consensus and variants
        all_consensus = allConsensus(runArtic.out[0].collect())
        tmp = runArtic.out.pass_vcf.toList().transpose().toList() // surely theres another way?
        all_variants = allVariants(tmp, reference)
        // genotype summary
        if (params.genotype_variants) {
            genotype_summary = genotypeSummary(
                runArtic.out.merged_gvcf, runArtic.out.bam, ref_variants).collect()
            combined_genotype_summary = combineGenotypeSummaries(genotype_summary)
        } else {
            genotype_summary = Channel.fromPath("$projectDir/data/OPTIONAL_FILE")
        }
        // nextclade
        clades = nextclade(all_consensus[0], reference, primers)
        // pangolin
        lineages = pangolin(all_consensus[0])
        // report
        html_doc = report(
            runArtic.out.depth_stats.collect(),
            read_summaries.collect(), 
            clades.collect(),
            lineages.collect(),
            genotype_summary.collect(),
            runArtic.out.vcf_stats.collect(),
            all_consensus[1])
        results = all_consensus[0].concat(all_consensus[1], all_variants[0].flatten(), 
            html_doc, combined_genotype_summary)
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

    if (!valid_scheme_versions.any { it == params.scheme_version}) {
        println("`--scheme_version` should be one of: $valid_scheme_versions, for `--scheme_name`: $params.scheme_name")
        exit 1
    }

    if (params.scheme_name == "spike-seq" && !params.genotype_variants) {
        println("`--genotype_variants` is required for scheme: 'spike-seq'")
        exit 1
    }

    if (params.samples && params.detect_samples) {
        println("Select either `--samples` or `--detect_samples`, not both")
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
        params._max_len = params.max_len
        params.remove('max_len')
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
    
    // check genotype variants
    if (params.genotype_variants) {
        if (params.genotype_variants == true) {
            ref_variants = file(
                "${scheme_directory}/${params.full_scheme_name}/${params.scheme_name}.vcf",
                type:'file', checkIfExists:true)
        } else {
            ref_variants = file(params.genotype_variants, type:'file', checkIfExists:true)
        }
    } else {
        ref_variants = Channel.fromPath("$projectDir/data/OPTIONAL_FILE")
    }

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
            // check if we need to autodetect sample names
            if (params.detect_samples) {
                sample_sheet = Channel
                    .fromPath(valid_barcode_dirs)
                    .map { path -> tuple(
                        path.baseName, 
                        file("${path}{**,.}/*.{fastq,fastq.gz,fq,fq.gz}")) }
                    .filter { pathTuple -> !pathTuple[1].isEmpty() }
                    .map { pathTuple -> tuple(
                        pathTuple[0],
                        sampleFromHeader(pathTuple[1][0])) }
            } else {
                // just map directory name to self
                sample_sheet = Channel
                    .fromPath(valid_barcode_dirs)
                    .filter(~/.*barcode[0-9]{1,3}$/)  // up to 192
                    .map { path -> tuple(path.baseName, path.baseName) }
            }
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
    results = pipeline(samples, scheme_directory, reference, primers, ref_variants)
    output(results)
}
