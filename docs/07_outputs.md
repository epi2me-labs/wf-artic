Output files may be aggregated including information for all samples or provided per sample. Per-sample files will be prefixed with respective aliases and represented below as {{ alias }}.

| Title | File path | Description | Per sample or aggregated |
|-------|-----------|-------------|--------------------------|
| Workflow report | ./wf-artic-report.html | Report for all samples. | aggregated |
| Consensus sequences | ./all_consensus.fasta | Final consensus sequences for all samples in the analysis. | aggregated |
| Pangolin results | ./lineage_report.csv | Pangolin results for each of the samples in the analysis. | aggregated |
| Nextclade results | ./nextclade.json | Nextclade results for each of the samples in the analysis. | aggregated |
| Coverage data | ./all_depth.txt | Coverage of the reference genome in 20 base windows in all the samples in the analysis. | aggregated |
| Variants | ./{{ alias }}.pass.named.vcf.gz | A VCF file containing high confidence variants in the sample when compared to the reference. | per-sample |
| Variants index | ./{{ alias }}.pass.named.vcf.gz.tbi | An index file for the variants. | per-sample |
| Alignments | ./{{ alias }}.primertrimmed.rg.sorted.bam | A BAM file containing the reads for the sample aligned to the reference. | per-sample |
| Alignments index | ./{{ alias }}.primertrimmed.rg.sorted.bam.bai | An index file for the alignments. | per-sample |
