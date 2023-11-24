The wf-artic workflow implements a slightly modified ARTIC FieldBioinformatics
workflow for the purpose of preparing consensus sequences from SARS-CoV-2
genomes that have been DNA sequenced using a pooled tiling amplicon strategy.

The workflow consumes a folder containing demultiplexed sequence reads as
prepared by either MinKNOW or Guppy. The workflow needs to know the primer
scheme that has been used during genome amplification and library preparation
e.g. `ARTIC/V3` or `ONT_Midnight/V1`. Other parameters can be specified too e.g.
assign sample names to the barcodes or to adjust the length distribution of
acceptable amplicon sequences. The Medaka variant model is selected based on the
provided basecaller configuration (using the parameter `--basecaller_cfg`), or
alternatively the Medaka model can be provided directly via the `--medaka_variant_model`
parameter.