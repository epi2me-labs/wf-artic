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
    log.info '    --samples        FILE    (Optional) tab-separated file with columns `barcode`'
    log.info '                             and `sample_name` indicating correspondence between '
    log.info '                             barcodes and sample names.'
    log.info ''

    return
}


process maybeDemultiplex {
    label "artic"
    cpus 1
    input:
        file directory
    output:
        file "out_directory"
    """
    #!/usr/bin/env python

    import os
    import pysam

    found = False
    for d in os.listdir("$directory"):
        if d.startswith("barcode"):
            for f in os.listdir(os.path.join("$directory", d)):
                fpath = os.path.join(d, f)
                try:
                    read = next(pysam.FastxFile(fpath))
                except Exception as e:
                    print(e)
                    print("{} was not readable as fastq".format(fpath))
                    pass
                else:
                    print("{} looks like a fastq".format(fpath))
                    found = True
                    break
        if found:
            break
    if found:
        os.symlink("$directory", "out_directory")
    else:
        # TODO: run demultiplexing
        raise RuntimeError("No 'barcode' directories found with valid .fastq")
    """
}


process preArticQC {
    label "artic"
    cpus 1
    input:
        file directory
    output:
        file "out_directory"
    """
    echo $directory
    ln -s $directory out_directory
    """
}


process runArtic {

    label "artic"
    cpus 2
    input:
        file directory
    output:
        file "*.consensus.fasta"

    """
    echo "=========="
    echo $directory
    echo "=========="
    artic guppyplex --skip-quality-check \
        --min-length $params.min_len --max-length $params.max_len \
        --directory $directory --prefix $params.prefix \
        && echo " - artic guppyplex finished"
    # the output of the above will be...
    READFILE=$params.prefix"_"$directory".fastq"

    artic minion --medaka --normalise 200 --threads $task.cpus \
        --read-file \$READFILE $params.scheme \
        --medaka-model $params.medaka_model \
        $directory \
        && echo " - artic minion finished"
    # output of the above will be...
    sed -i "s/^>.*/>$directory/" $directory".consensus.fasta"
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
        file "bc_cons_*.fastq"
        //file (optional) barcode_to_sample_key
    output:
        file "all_consensus.fasta"
    """
    cat bc_cons_*.fastq > all_consensus.fasta
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
        source_directory
    main:
        demultiplexed = maybeDemultiplex(source_directory)
        barcode_dirs = channel.fromPath(params.fastq + "/barcode*/", type: 'dir') 
        barcode_dirs.view { "value: $it" }
        qc_results = preArticQC(barcode_dirs)
        artic_consensus = runArtic(barcode_dirs)
        consensus = allConsensus(artic_consensus)
    emit:
        consensus
}

// entrypoint workflow
workflow {
    barcode_dirs = channel.fromPath(params.fastq, type: 'dir')
    barcode_dirs.view { "value: $it" }
    results = pipeline(barcode_dirs)
    output(results)
}
