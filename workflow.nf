#!/usr/bin/env extflow

nextflow.enable.dsl = 2

params.help = ""

if(params.help) {
    log.info ''
    log.info 'Workflow template'
    log.info ''
    log.info 'Usage: '
    log.info '    nextflow run workflow.nf [options]'
    log.info ''
    log.info 'Script Options: '
    log.info '    --fastq        DIR     Path to FASTQ directory'
    log.info '    --out_dir      DIR     Path for output'
    log.info '    --samples      FILE    (Optional) tab-separated file with columns `barcode`'
    log.info '                           and `sample_name` indicating correspondence between '
    log.info '                           barcodes and sample names.'
    log.info ''

    return
}


process maybeDemultiplex {
    // Check the input path, demultiplex if necessary, emit barcode folders
    label "artic"
    cpus 1    
    input:
        file reads
    output:
        stdout

    """
    # TODO

    """
}

process QC {
    label "artic"
    cpus 1
    input:
        file directory
    output:
        file out_directory
    """
    # TODO: anything?
    """
}


process runArtic {

    label "artic"
    cpus 2
    input:
        file directory
    output:
        file consensus.fasta

    """
    # TODO - straighten this out
    artic guppyplex --skip-quality-check \
        --min-length $min_len --max-length $max_len \
        --directory $directory --prefix $prefix >>$log_file 2>&1 \
        && echo " - artic guppyplex finished"

    artic minion --medaka --normalise 200 --threads $threads \
        --scheme-directory $working_dir/artic-ncov2019/primer_schemes \
        --read-file $read_file $scheme \
        --medaka-model $medaka_model \
        $prefix >>$log_file 2>&1 \
        && echo " - artic minion finished"

    # TODO: copied from notebook, calculates depth per amplicon pool
    for i in (1,2):
        input_bam = "{}.primertrimmed.rg.sorted.bam".format(prefix)
        tag = "nCoV-2019_{}".format(i)
        bam = "{}.{}".format(input_bam, tag)
        if os.path.isfile(input_bam):
            !samtools view -r $tag -b $input_bam > $bam
            !samtools index $bam
            !stats_from_bam $bam > $bam".stats"
            !coverage_from_bam -s 50 -p $bam $bam
            !rm -rf $bam $bam.bai
        else:
            !echo error " - No Artic output bam file found for primer set $i, Artic failed."
    """
}

process report {
    label "artic"
    cpus 1
    input:
        gather_some_stuff
    output:
        file report_file
    """
    #TODO make a report from all barcodes
    """
}

process all_consensus {
    label "artic"
    cpus 1
    input:
        gather_some_stuff
        file (optional) barcode_to_sample_key
    output:
        file all_consensus.fasta
    """
    # TODO: gather all consensus, concatenate whilst (optionally) renaming
    """


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
        reads
    main:
        seqs = maybeDemultiplex(reads).splitText()
    emit:
        seqs
}

// entrypoint workflow
workflow {
    reads = channel.fromPath(params.reads, checkIfExists:true)
    results = pipeline(reads)
    output(results)
}
