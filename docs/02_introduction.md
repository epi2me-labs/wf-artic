> This workflow is deprecated. It is not longer proactively 
> maintained; and reactive maintainance is sporadic and low priority for our
> developers. Although the repository is still available, it is now unsupported
> and we do not recommend its use. Please contact support@nanoporetech.com for
> help with your application.
>
> For a more actively maintained workflow, please consider using the upstream
> [fieldbioinformatics](https://github.com/artic-network/fieldbioinformatics/)
> project directly.

The wf-artic workflow implements a slightly modified ARTIC FieldBioinformatics
workflow for the purpose of preparing consensus sequences from SARS-CoV-2
genomes that have been DNA sequenced using a pooled tiling amplicon strategy.

The workflow consumes a folder containing demultiplexed sequence reads as
prepared by either MinKNOW or Guppy. The workflow needs to know the primer
scheme that has been used during genome amplification and library preparation
e.g. `ARTIC/V3` or `ONT_Midnight/V1`. Other parameters can be specified too e.g.
assign sample names to the barcodes or to adjust the length distribution of
acceptable amplicon sequences.
