### Input Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| fastq | string | FASTQ files to use in the analysis. | This accepts one of three cases: (i) the path to a single FASTQ file; (ii) the path to a top-level directory containing FASTQ files; (iii) the path to a directory containing one level of sub-directories which in turn contain FASTQ files. In the first and second case, a sample name can be supplied with `--sample`. In the last case, the data is assumed to be multiplexed with the names of the sub-directories as barcodes. In this case, a sample sheet can be provided with `--sample_sheet`. |  |
| analyse_unclassified | boolean | Analyse unclassified reads from input directory. By default the workflow will not process reads in the unclassified directory. | If selected and if the input is a multiplex directory the workflow will also process the unclassified directory. | False |


### Primer Scheme Selection

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| scheme_name | string | Primer scheme name. | This should be set to `SARS-CoV-2`, or `spike-seq` or your custom scheme name. This affects the choice of scheme versions you can use. The only scheme versions compatible with `spike-seq` are `ONT/V1` and `ONT/V4.1` | SARS-CoV-2 |
| scheme_version | string | Primer scheme version. | This is the version of the primer scheme to use, more details about primer shemes can be found [here](https://labs.epi2me.io/ont-midnight-scheme-update/). | ARTIC/V3 |
| custom_scheme | string | Path to a custom scheme. | If you have a custom primer scheme you can enter the details here. This must be the full path to the directory containing your appropriately named scheme bed and fasta files; <SCHEME_NAME>.bed and <SCHEME_NAME>.fasta. More details [here](https://labs.epi2me.io/ont-midnight-scheme-update/). |  |


### Sample Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| sample_sheet | string | A CSV file used to map barcodes to sample aliases. The sample sheet can be provided when the input data is a directory containing sub-directories with FASTQ files. | The sample sheet is a CSV file with, minimally, columns named `barcode` and `alias`. Extra columns are allowed. A `type` column is required for certain workflows and should have the following values; `test_sample`, `positive_control`, `negative_control`, `no_template_control`. |  |
| sample | string | A single sample name for non-multiplexed data. Permissible if passing a single .fastq(.gz) file or directory of .fastq(.gz) files. |  |  |


### Output Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| out_dir | string | Directory for output of all workflow results. |  | output |


### Reporting Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| report_depth | integer | Min. depth for percentage coverage. (e.g. 89% genome covered at > `report_depth`) |  | 100 |
| report_clade | boolean | Show results of Nextclade analysis in report. |  | True |
| report_coverage | boolean | Show genome coverage traces in report. |  | True |
| report_lineage | boolean | Show results of Pangolin analysis in report. |  | True |
| report_variant_summary | boolean | Show variant information in report. |  | True |


### Advanced Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| artic_threads | number | Number of CPU threads to use per artic task. | The total CPU resource used by the workflow is constrained by the executor configuration. | 4 |
| pangolin_threads | number | Number of CPU threads to use per pangolin task. | The total CPU resource used by the workflow is constrained by the executor configuration. | 4 |
| genotype_variants | string | Report genotyping information for scheme's known variants of interest, optionally provide file path as argument. |  |  |
| list_schemes | boolean | List primer schemes and exit without running analysis. |  | False |
| min_len | number | Minimum read length (default: set by scheme). |  |  |
| max_len | number | Maximum read length (default: set by scheme). |  |  |
| max_softclip_length | integer | Remove reads with alignments showing large soft clipping |  |  |
| update_data | boolean | Update Pangolin and Nextclade data at runtime. |  | True |
| pangolin_options | string | Pass options to Pangolin, for example "--analysis-mode fast --min-length 26000". |  |  |
| nextclade_data_tag | string | The tag of the nextclade data packet |  |  |
| normalise | integer | Depth ceiling for depth of coverage normalisation |  | 200 |
| override_basecaller_cfg | string | Override auto-detected basecaller model that processed the signal data; used to select an appropriate Medaka model. | Per default, the workflow tries to determine the basecall model from the input data. This parameter can be used to override the detected value (or to provide a model name if none was found in the inputs). However, users should only do this if they know for certain which model was used as selecting the wrong option might give sub-optimal results. A list of recent models can be found here: https://github.com/nanoporetech/dorado#DNA-models. |  |


### Miscellaneous Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| lab_id | string | Laboratory identifier, used in reporting. |  |  |
| testkit | string | Test kit identifier, used in reporting. |  |  |


