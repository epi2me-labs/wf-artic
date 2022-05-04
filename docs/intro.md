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













