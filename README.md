# ARTIC SARS-CoV-2 Workflow

This repository contains a [nextflow](https://www.nextflow.io/) workflow for
running the ARTIC SARS-CoV-2 workflow on multiplexed MinION, GridION, and
PromethION runs.
## Introduction

The wf-artic workflow implements a slightly modified ARTIC FieldBioinformatics
workflow for the purpose of preparing consensus sequences from SARS-CoV-2
genomes that have been DNA sequenced using a pooled tiling amplicon strategy.

The workflow consumes a folder containing demultiplexed sequence reads as
prepared by either MinKNOW or Guppy. The workflow needs to know the primer
scheme that has been used during genome amplication and library preparation
e.g. ARTIC/V3 or ONT_Midnight/V1. Other parameters can be specified to e.g.
assign sample names to the barcodes or to adjust the length distribution of
acceptable amplicon sequences.

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













## Quickstart

The workflow uses [nextflow](https://www.nextflow.io/) to manage compute and
software resources, as such nextflow will need to be installed before attempting
to run the workflow.

The workflow can currently be run using either
[Docker](https://www.docker.com/products/docker-desktop) or
[conda](https://docs.conda.io/en/latest/miniconda.html) to provide isolation of
the required software. Both methods are automated out-of-the-box provided
either docker of conda is installed.

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

* [medaka](https://www.github.com/nanoporetech/medaka)
* [artic](https://github.com/artic-network/fieldbioinformatics)
* [nextflow](https://www.nextflow.io/)
* [docker](https://www.docker.com/products/docker-desktop)
* [conda](https://docs.conda.io/en/latest/miniconda.html)
