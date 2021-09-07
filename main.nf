#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { fastq_ingress } from './lib/fastqingress' 

valid_schemes = ["SARS-CoV-2", "spike-seq"]
valid_scheme_versions = ["V1", "V2", "V3", "V4", "V4.1", "V1200"]

if (params.scheme_name == "spike-seq") {
    valid_scheme_versions = ["V1", "V4.1"]
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
    --max_softclip_length       INT     Maximum alignment overhang length, to remove possibly chimeric reads (default: 0)
    --scheme_name               STR     Scheme to use ($valid_schemes) (default: SARS-CoV-2)
    --scheme_version            STR     Primer scheme version ($valid_scheme_versions)
                                        (default: $params.scheme_version)
    --report_depth              INT     Min. depth for percentage coverage (default: $params.report_depth)
                                        (e.g. 89% genome covered at > `report_depth`)
                                        indicating correspondence between
    --genotype_variants         FILE    Report genotyping information for scheme's known variants of interest,
                                        optionally provide file path as argument.
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
    fastcat -s ${sample_name} -r ${sample_name}.stats -x ${directory} > /dev/null
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
        tuple path("${sample_name}.primertrimmed.rg.sorted.bam"), path("${sample_name}.primertrimmed.rg.sorted.bam.bai"), emit: bam
    """
    run_artic.sh \
        ${sample_name} ${directory} ${params._min_len} ${params._max_len} \
        ${params.medaka_model} ${params.full_scheme_name} \
        ${task.cpus} ${params._max_softclip_length}
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


process get_versions {
    label "artic"
    cpus 1
    output:
        path "versions.txt"
    script:
    """
    medaka --version | sed 's/ /,/' >> versions.txt
    minimap2 --version | sed 's/^/minimap2,/' >> versions.txt
    bcftools --version | head -n 1 | sed 's/ /,/' >> versions.txt
    samtools --version | head -n 1 | sed 's/ /,/' >> versions.txt
    nextclade --version | sed 's/^/nextclade,/' >> versions.txt
    artic --version | sed 's/ /,/' >> versions.txt
    """
}


process report {
    label "artic"
    cpus 1
    input:
        path "depth_stats/*"
        path "read_stats/*"
        path "nextclade.json"
        path "pangolin.csv"
        path "genotypes/*"
        path "vcf_stats/*"
        path "consensus_status.txt"
        path "versions/*"
    output:
        path "wf-artic-report.html"
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
        --bcftools_stats vcf_stats/* $genotype \
        --versions versions
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
        file "nextclade_dataset"
    output:
        file "nextclade.json"
    """
    cp -L reference.fasta ref.fasta
    scheme_to_nextclade.py $scheme_bed ref.fasta primers.csv
    nextclade run \
        --input-fasta consensus.fasta \
        --reference nextclade_dataset/reference.fasta \
        --input-pcr-primers primers.csv \
        --input-tree nextclade_dataset/tree.json \
        --input-qc-config nextclade_dataset/qc.json \
        --input-gene-map nextclade_dataset/genemap.gff \
        --output-json nextclade.json \
        --jobs 1
    """
}


process pangolin {
    label "pangolin"
    cpus 1
    input:
        path "consensus.fasta"
    output:
        path "lineage_report.csv", emit: report
        path "pangolin.version", emit: version
    """
    pangolin --version 2>&1 | sed 's/ /,/' > pangolin.version
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
        nextclade_dataset
    main:
        software_versions = get_versions()
        combined_genotype_summary = Channel.empty()
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
        clades = nextclade(
            all_consensus[0], reference, primers, nextclade_dataset)
        // pangolin
        pangolin(all_consensus[0])
        software_versions = software_versions.mix(pangolin.out.version)
        // report
        html_doc = report(
            runArtic.out.depth_stats.collect(),
            read_summaries.collect(), 
            clades.collect(),
            pangolin.out.report.collect(),
            genotype_summary.collect(),
            runArtic.out.vcf_stats.collect(),
            all_consensus[1],
            software_versions.collect())
        results = all_consensus[0].concat(all_consensus[1], all_variants[0].flatten(),
            runArtic.out.bam.flatten(), html_doc, combined_genotype_summary)
        results.view()
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
    if (!params.max_softclip_length) {
        params.remove('max_softclip_length')
        params._max_softclip_length = 0
    }
    else{
        params._max_softclip_length = params.max_softclip_length
        params.remove('max_softclip_length')
    }
    println("")
    println("Parameter summary")
    println("=================")
    params.each { it -> println("    $it.key: $it.value") }
    println("")

    params.full_scheme_name = params.scheme_name + "/" + params.scheme_version

    schemes = projectDir + '/data/primer_schemes'
    scheme_directory = file(schemes, type: 'dir', checkIfExists:true)
    nextclade = projectDir + '/data/nextclade/'
    nextclade_dataset = file(nextclade, type: 'dir', checkIfExists:true)

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

    // check fastq dataset and run workflow
    samples = fastq_ingress(
        params.fastq, workDir, params.samples, params.sanitize_fastq)
    results = pipeline(samples, scheme_directory, reference, 
        primers, ref_variants, nextclade_dataset)
    output(results)
}
