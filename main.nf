#!/usr/bin/env nextflow

import groovy.json.JsonBuilder
nextflow.enable.dsl = 2

include { fastq_ingress } from './lib/fastqingress'


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
        tuple file(directory), val(sample_id), val(type)
    output:
        file "${sample_id}.stats"
    """
    fastcat -s ${sample_id} -r ${sample_id}.stats -x ${directory} > /dev/null
    """
}


process runArtic {
    label "artic"
    cpus 2
    input:
        tuple file(directory), val(sample_id), val(type)
        file "primer_schemes"
    output:
        path "${sample_id}.consensus.fasta", emit: consensus
        path "${sample_id}.depth.txt", emit: depth_stats
        path "${sample_id}.pass.named.stats", emit: vcf_stats
        tuple(
            val(sample_id),
            path("${sample_id}.pass.named.vcf.gz"),
            path("${sample_id}.pass.named.vcf.gz.tbi"),
            emit: pass_vcf)
        tuple(
            val(sample_id),
            path("${sample_id}.merged.gvcf.named.vcf.gz"),
            path("${sample_id}.merged.gvcf.named.vcf.gz.tbi"),
            emit: merged_gvcf)
        tuple(
            val(sample_id),
            path("${sample_id}.primertrimmed.rg.sorted.bam"),
            path("${sample_id}.primertrimmed.rg.sorted.bam.bai"),
            emit: primertrimmed_bam)
        tuple(
            val(sample_id),
            path("${sample_id}.trimmed.rg.sorted.bam"),
            path("${sample_id}.trimmed.rg.sorted.bam.bai"),
            emit: trimmed_bam)
    """
    run_artic.sh \
        ${sample_id} ${directory} ${params._min_len} ${params._max_len} \
        ${params.medaka_model} ${params.full_scheme_name} \
        ${task.cpus} ${params._max_softclip_length}
    bcftools stats ${sample_id}.pass.named.vcf.gz > ${sample_id}.pass.named.stats 
    """
}


process genotypeSummary {
    // Produce a genotype summary spreadsheet
    label "artic"
    cpus 1
    input:
        tuple val(sample_id), file(vcf), file(tbi), file(bam), file(bam_index)
        file "reference.vcf"
    output:
        file "*genotype.csv"
    script:
        def lab_id = params.lab_id ? "--lab_id ${params.lab_id}" : ""
        def testkit = params.testkit ? "--testkit ${params.testkit}" : ""
    """
    genotype_summary.py \
        -b $bam \
        -v $vcf \
        -d reference.vcf \
        --sample $sample_id \
        $lab_id \
        $testkit \
        -o ${csvName}.genotype.csv
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


process getVersions {
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


process getParams {
    label "wfplasmid"
    cpus 1
    output:
        path "params.json"
    script:
        def paramsJSON = new JsonBuilder(params).toPrettyString()
    """
    # Output nextflow params object to JSON
    echo '$paramsJSON' > params.json
    """
}


process telemetry {
    label "artic"
    cpus 1
    input:
        tuple val(sample_id), file(bams), file(bais), file(vcfs), file(tbis)
        path scheme_bed
        path reference
    output:
        path "telemetry.json", emit: json
    script:
        def samples = sample_id.join(' ')
    """
    output_telemetry.py \
        telemetry.json \
        --scheme_name $params.scheme_name \
        --scheme_bed $scheme_bed \
        --reference $reference \
        --samples $samples \
        --alignments $bams \
        --calls $vcfs
    """
}


process report {
    label "artic"
    cpus 1
    input:
        path "depth_stats/*"
        path "read_stats/*"
        path "nextclade.json"
        path nextclade_errors
        path "pangolin.csv"
        path "genotypes/*"
        path "vcf_stats/*"
        path "consensus_status.txt"
        path "versions/*"
        path "params.json"
        path "consensus_fasta"
        path "telemetry.json"
        val samples
        val types
    output:
        path "wf-artic-*.html"
        path "*.json"
    script:
    // when genotype_variants is false the channel contains a mock file
    def report_name = "wf-artic-" + params.report_name + '.html'
    def genotype = params.genotype_variants ? "--genotypes genotypes/*" : ""
    def nextclade = params.report_clade as Boolean ? "--nextclade nextclade.json" : ""
    def pangolin = params.report_lineage as Boolean ? "--pangolin pangolin.csv" : ""
    def coverage = params.report_coverage as Boolean ? "" : "--hide_coverage"
    def var_summary = params.report_variant_summary as Boolean ? "" : "--hide_variants"
    def debug = params.report_detailed as Boolean ? "--telemetry telemetry.json" : "--hide_debug"
    """
    echo "$pangolin"
    echo "$nextclade"
    report.py \
        consensus_status.txt $report_name \
        $pangolin $coverage $var_summary \
        $nextclade $debug \
        --nextclade_errors $nextclade_errors \
        --revision $workflow.revision \
        --commit $workflow.commitId \
        --min_len $params._min_len \
        --max_len $params._max_len \
        --report_depth $params.report_depth \
        --depths depth_stats/* \
        --summaries read_stats/* \
        --bcftools_stats vcf_stats/* $genotype \
        --versions versions \
        --params params.json \
        --consensus_fasta $consensus_fasta \
        --samples $samples \
        --types $types
    """
}


process report_no_data {
    label "artic"
    cpus 1
    input:
        path "versions/*"
        val error
        path "params.json"
    output:
        path "wf-artic-*.html"
        path "*.json", optional: true
    script:
    // when genotype_variants is false the channel contains a mock file
    def report_name = "wf-artic-" + params.report_name + '.html'
    def error_message = error
    """
    report_error.py \
        --output $report_name \
        --revision $workflow.revision --params params.json --commit $workflow.commitId \
        --versions versions --error_message \"$error_message\"
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
        tuple val(sample_id), file(vcfs), file(tbis)
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
        file "*.errors.csv"
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
    pangolin --all-versions 2>&1 | sed 's/: /,/' > pangolin.version
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
        software_versions = getVersions()
        workflow_params = getParams()
        combined_genotype_summary = Channel.empty()
        scheme_directory = copySchemeDir(scheme_directory)
        if ((samples.getClass() == String) && (samples.startsWith("Error"))){
            samples = channel.of(samples)
            html_doc = report_no_data(
                software_versions.collect(),
                samples,
                workflow_params)
            results = html_doc[0].concat(html_doc[1])
        } else {
            read_summaries = preArticQC(samples)
            artic = runArtic(samples, scheme_directory)
            // collate consensus and variants
            artic.consensus.view()
            all_consensus = allConsensus(artic.consensus.collect())
            all_variants = allVariants(
                artic.pass_vcf.toList().transpose().toList(), reference)
            // genotype summary
            if (params.genotype_variants) {
                genotype_summary = genotypeSummary(
                    artic.merged_gvcf.join(artic.primertrimmed_bam), ref_variants)
                combined_genotype_summary = combineGenotypeSummaries(
                    genotype_summary.collect())
            } else {
                genotype_summary = Channel.fromPath("$projectDir/data/OPTIONAL_FILE")
            }
            // nextclade
            clades = nextclade(
                all_consensus[0], reference, primers, nextclade_dataset)
            // pangolin
            pangolin(all_consensus[0])
            software_versions = software_versions.mix(pangolin.out.version)
            // telemetry
            telemetry_output = telemetry(
                artic.trimmed_bam.join(artic.pass_vcf).toList().transpose().toList(),
                primers,
                reference)
            // report
            
            html_doc = report(
                artic.depth_stats.collect(),
                read_summaries.collect(), 
                clades[0].collect(),
                clades[1].collect(),
                pangolin.out.report.collect(),
                genotype_summary.collect(),
                artic.vcf_stats.collect(),
                all_consensus[1],
                software_versions.collect(),
                workflow_params,
                all_consensus[0],
                telemetry_output,
                // sample_ids
                samples.map{ it -> it[1]}.toList().map{ it.join(' ')},
                // sample types
                samples.map{ it -> it[2]}.toList().map{ it.join(' ')}
                )
            results = all_consensus[0].concat(
                telemetry.out.json,
                all_consensus[1],
                all_variants[0].flatten(),
                clades[0],
                artic.primertrimmed_bam.flatMap { it -> [ it[1], it[2] ] },
                html_doc[0],
                html_doc[1],
                combined_genotype_summary,
                pangolin.out.report)
            }
    emit:
        results
}


// entrypoint workflow
WorkflowMain.initialise(workflow, params, log)

valid_schemes = ["SARS-CoV-2", "spike-seq"]
valid_scheme_versions = ["V1", "V2", "V3", "V4", "V4.1", "V1200"]

if (params.scheme_name == "spike-seq") {
    valid_scheme_versions = ["V1", "V4.1"]
}

workflow {

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
