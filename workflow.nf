#!/usr/bin/env extflow

nextflow.enable.dsl = 2

params.help = ""
params.out_dir = "output"
params.prefix = "artic"
params.min_len = "300"
params.max_len = "700"
params.medaka_model = "r941_min_high_g360"
params.scheme = "SARS-CoV-2/V3"



if(params.help) {
    log.info ''
    log.info 'Workflow template'
    log.info ''
    log.info 'Usage: '
    log.info '    nextflow run workflow.nf [options]'
    log.info ''
    log.info 'Script Options: '
    log.info '    --fastq          DIR     Path to FASTQ directory'
    log.info '    --out_dir        DIR     Path for output'
    log.info '    --prefix         STR     Output filename prefix'
    log.info '    --medaka_model   STR     Medaka model name'
    log.info '    --min_length     INT     Minimum read length'
    log.info '    --max_length     INT     Maximum read length'
    log.info '    --samples        FILE    CSV file with columns named `barcode` and `sample_name`' 
    log.info '                             indicating correspondence between'
    log.info '                             barcodes and sample names.'
    log.info ''

    return
}


process preArticQC {
    label "artic"
    cpus 1
    input:
        tuple file(directory), val(sample_name) 
    output:
        tuple file("out_directory"), val(sample_name) 
    """
    echo $sample_name
    ln -s $directory out_directory
    """
}


process runArtic {

    label "artic"
    cpus 2
    input:
        tuple file(directory), val(sample_name)
    output:
        file("${sample_name}.consensus.fasta")
        tuple file("${sample_name}.pass.named.vcf.gz"), file("${sample_name}.pass.named.vcf.gz.tbi")

    """
    echo "=========="
    echo $directory, $sample_name
    echo "=========="
    # name everything by the sample rather than barcode
    ln -s $directory $sample_name
    artic guppyplex --skip-quality-check \
        --min-length $params.min_len --max-length $params.max_len \
        --directory $sample_name --prefix $sample_name \
        && echo " - artic guppyplex finished"
    # the output of the above will be...
    READFILE=$sample_name"_"$sample_name".fastq"

    artic minion --medaka --normalise 200 --threads $task.cpus \
        --read-file \$READFILE $params.scheme \
        --medaka-model $params.medaka_model \
        $sample_name \
        && echo " - artic minion finished"
    # add the sample name to the pass VCF
    zcat $sample_name".pass.vcf.gz" | sed "s/SAMPLE/$sample_name/" | bgzip > $sample_name".pass.named.vcf.gz"
    bcftools index -t $sample_name".pass.named.vcf.gz" 
    # output of the above will be...
    sed -i "s/^>.*/>$sample_name/" $sample_name".consensus.fasta"
    """
}


process postArticQC {
    label "artic"
    cpus 1
    """
    echo "Nothing to see here"
    """
}

//process report {
//    label "artic"
//    cpus 1
//    input:
//        gather_some_stuff
//    output:
//        file report_file
//    """
//    #TODO make a report from all barcodes
//    """
//}


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
    main:
        preArticQC(samples)
        runArtic(samples)
        all_consensus = allConsensus(runArtic.out[0].collect())
        // surely theres another way?
        tmp = runArtic.out[1].toList().transpose().toList()
        all_variants = allVariants(tmp)
        results = all_consensus.concat(all_variants)
    emit:
        results
}

// entrypoint workflow
workflow {

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
    samples.view()
    results = pipeline(samples)
    output(results)
}
