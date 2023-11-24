The pipeline is largely a wrapper around the [Artic Network](https://artic.network/) [Field Bioinformatics](https://github.com/artic-network/fieldbioinformatics) analysis package.

### 1. Concatenates input files and generate per read stats.

The [fastcat/bamstats](https://github.com/epi2me-labs/fastcat) tool is used to concatenate multifile samples to be processed by the workflow. It will also output per read stats including average read lengths and qualities. Reads are additionally filtered for sequence length and quality characteristics.

### 2. Mapping and primer trimming (Artic)

Concatenated reads are mapped to the reference SARS-CoV-2 genome using [minimap2](https://github.com/lh3/minimap2). A primer scheme-specific BED file is used to identify the regions of
the mapped sequences that correspond to synthetic sequences (primers) - these regions are clipped to ensure that sequences are entirely of biological origin.

### 3. Variant calling and consensus generation (Artic)

The retained sequences are used to prepare a consensus sequence that is then polished using [Medaka](https://github.com/nanoporetech/medaka) and variant calling is performed to produce a VCF file of genetic differences relative to the reference genome.

### 4. Lineage/clade assignment

The consensus sequence is annotated for virus clade information using [NextClade](https://clades.nextstrain.org/), and strain assignment is performed using [Pangolin](https://github.com/cov-lineages/pangolin).