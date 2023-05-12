# ARTIC SARS-CoV-2 Workflow

This repository contains a [Nextflow](https://www.nextflow.io/) workflow for
running the ARTIC SARS-CoV-2 workflow on multiplexed MinION, GridION, and
PromethION runs.




## Introduction

The wf-artic workflow implements a slightly modified ARTIC FieldBioinformatics
workflow for the purpose of preparing consensus sequences from SARS-CoV-2
genomes that have been DNA sequenced using a pooled tiling amplicon strategy.

The workflow consumes a folder containing demultiplexed sequence reads as
prepared by either MinKNOW or Guppy. The workflow needs to know the primer
scheme that has been used during genome amplification and library preparation
e.g. ARTIC/V3 or ONT_Midnight/V1. Other parameters can be specified too e.g.
assign sample names to the barcodes or to adjust the length distribution of
acceptable amplicon sequences. The Medaka variant model is selected based on the
provided basecaller configuration (using the parameter `--basecaller_cfg`), or
alternatively the Medaka model can be provided directly via the `--medaka_variant_model`
parameter.

DNA sequences in FASTQ format are aggregated, filtered for sequence length and
quality characteristics and are mapped to the reference SARS-CoV-2 genome using
minimap2. A primer-scheme specific bed file is used to identify the regions of
the mapped sequences that correspond to synthetic sequences (primers) - these
regions are clipped to ensure that sequences are entirely of biological origin.
The retained sequences are used to prepare a consensus sequence that is then
polished using Medaka and variant calling is performed to produce a VCF file of
genetic differences relative to the reference genome. The consensus sequence is
annotated for virus clade information using NextClade and a strain assignment
is performed using Pangolin.

The completed analysis is summarised in an HTML format report that summarises
key information that includes number of reads, depth of coverage information
per amplicon and both the Nextclade and Pangolin information.

More information can be found in these two blog posts:
* [SARS-CoV-2 Midnight Analysis](https://labs.epi2me.io/sarscov2-midnight-analysis/)
* [Midnight Scheme Update](https://labs.epi2me.io/ont-midnight-scheme-update/)




## Quickstart

The workflow uses [Nextflow](https://www.nextflow.io/) to manage compute and
software resources, as such Nextflow will need to be installed before attempting
to run the workflow.

The workflow can currently be run using either
[Docker](https://www.docker.com/products/docker-desktop) or
[Singularity](https://docs.sylabs.io/guides/latest/user-guide/) to provide isolation of
the required software. Both methods are automated out-of-the-box provided
either Docker or Singularity is installed.

It is not required to clone or download the git repository in order to run the workflow.
For more information on running EPI2ME Labs workflows [visit our website](https://labs.epi2me.io/wfindex).

**Workflow options**

To obtain the workflow, having installed `nextflow`, users can run:

```
nextflow run epi2me-labs/wf-artic --help
```

to see the options for the workflow.

**Workflow outputs**

The primary outputs of the workflow include:

* a [FASTA](https://en.wikipedia.org/wiki/FASTA) file containing the consensus sequence for all samples,
* a [VCF](https://en.wikipedia.org/wiki/Variant_Call_Format) file sample all samples,
* an HTML report document detailing QC metrics and the primary findings of the workflow.




## Useful links

* [Medaka](https://www.github.com/nanoporetech/medaka)
* [Artic](https://github.com/artic-network/fieldbioinformatics)
* [Nextflow](https://www.nextflow.io/)
* [Docker](https://www.docker.com/products/docker-desktop)
* [Singularity](https://docs.sylabs.io/guides/latest/user-guide/)
