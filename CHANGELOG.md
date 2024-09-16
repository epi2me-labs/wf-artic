# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [v1.2.2]
### Changed
- Reconcile workflow with wf-template v5.2.6

## [v1.2.1]
### Changed
- Reconcile workflow with wf-template v5.2.5

## [v1.2.0]
### Added
- Logging to file for artic commands (`{barcode}.artic.log.txt`)
- `--override_basecaller_cfg` parameter for cases where automatic basecall model detection fails or users wish to override the automatic choice
### Removed
- Unused `nextclade_version` parameter
- The `--medaka_variant_model` parameter as the appropriate Medaka model is now automatically determined from the input data
- The now redundant `--basecaller_cfg` parameter as its value is now automatically detected from the input data on a per-sample basis
### Changed
- Reconciled workflow with wf-template v5.2.1
- Updated Medaka to v1.12.1

## [v1.1.0]
### Updated
- Nextclade updated to V3
### Removed
- Nextclade V2 datasets

## [v1.0.0]
### Changes
- Documentation
- Defined resource requirments

## [v0.3.33]
### Changes
- Updates for cloud readiness
- Docs updated
### Removed
- Default local executor CPU and RAM limits

## [v0.3.32]
### Fixed
- reporting of sequence summaries crashing with `TypeError`
- sed error when writing to tmp

## [v0.3.31]
### Changed
- `--update_data` default is now `true` to get latest lineage data from Pangolin and Nextclade
### Updated
- Pangolin and Nextclade containers update to latest data

## [v0.3.30]
### Added
- New primer schemes; ARTIC V5.3.2 and Midnight IDT V2
### Updated
- Nextclade bump to 2.14.0
- Pangolin bump to 4.3
- For alt primers we take the largest region (Updated ARTIC code merges a canonical primer site with an alt site, producing an interval that encompasses both)

## [v0.3.29]
### Changed
- Example command to use demo data.
- Updated GitHub issue templates to force capture of more information.
- Bumped minimum required Nextflow version to 22.10.8
- Enum choices are enumerated in the `--help` output.
- Enum choices are enumerated as part of the error message when a user has selected an invalid choice.
- Now uses hard-coded values instead of `params.threads` in `fastqingress`.

## [v0.3.28]
### Fixed
- Fixed an error which was causing only a single sample to be processed, instead of all samples in the input FASTQ directory

## [v0.3.27]
### Added
- Configuration for running demo data in AWS

## [v0.3.26]
### Changed
- `--basecall_cfg` is now used to determine a suitable Medaka model, alternatively provide the name of a model with `--medaka_variant_model` to override automatic selection
- Update documentation

### Fixed
- `--report_detailed` was causing an error in report generation, as it requires the telemetry JSON which the workflow is no longer using. This parameter has been removed (https://github.com/epi2me-labs/wf-artic/issues/79)

### Added
- Configuration for running demo data in AWS

## [v0.3.25]
### Updated
- Implemented new fastqingress
- aplanat bumped to v0.6.16 to fix nextclade sequence view
- Updated to Oxford Nanopore Technologies PLC. Public License

## [v0.3.24]
### Fixed
- Custom primer schemes

## [v0.3.23]
### Changed
- medaka bumped to 1.7.2 (np-artic bump to 1.3.0_19)

### Added
- medaka models for kit14 - R10.4.1

### Fixed
- scheme directory bug in labs/windows (https://github.com/epi2me-labs/wf-artic/issues/73)

## [v0.3.22]
### Changed
- workflow-glue python psuedo package
- pangolin to 4.2
- nextclade to 2.11.0

### Removed
- spike-seq primer schemes
- telemtry generation

## [v0.3.21]
### Fixed
- Updated sample sheet to expect a file

## [v0.3.20]
### Changed
- Updated description in manifest

## [v0.3.19]
### Changed
- Fastqingress metadata map
- Nextclade < v2.0.0 no longer supported
- Documentation
- pangolin to 4.1.3
- nextclade to 2.8.0

### Fixed
- Metadata map handling

## [v0.3.18]
### Changes
- pangolin to v4.1.1
- nextclade data updated

## [v0.3.17]
### Changes
- Args parser for fastqingress
- Set out_dir option type to ensure output is written to correct directory on Windows

## [v0.3.16]
### Added
- `--pangolin_options` command line arg, use with quotes i.e.: "--analysis-mode fast"
- NEB VarSkip V2b primer scheme added

### Changes
- Better help text on cli

## [v0.3.15]
### Fixed
- `nextflow_schema.json` default primer scheme fix

### Changes
- New docs format

## [v0.3.14]
### Fixed
- Issue with custom schemes

## [v0.3.13]
### Changes
- Updated primer schemes
- Removed scheme downloading from np-artic
- V1200 Midnight naming changed:
- v1 of the kit (MRT.001.10) is called Midnight-ONT/V1
- v2 of the kit (MRT.001.20) is called Midnight-ONT/V2
- v3 of the kit (MRT.001.30) is called Midnight-ONT/V3
- It should be noted that primer 21_RIGHT only has one record in the bed, to represent the maximum possible region covered by the right primers for this amplicon
- pangolin bumped to 4.0.5
- nextclade bumped to 1.11.0

### Added
- Ability to specify custom primers with --custom_scheme
  - --min_len and --max_len must be specified when using this option
- --list_schemes option that lists available schemes and quits the workflow

## [v0.3.12]
### Changes
- Updated medaka default model to `r941_min_hac_variant_g507`
- Nextclade bumped to 1.10.3
- Help text for sample sheets updated as suggested in https://github.com/epi2me-labs/wf-artic/issues/41
- Output a combined depth file for later analysis

## [v0.3.11]
### Changes
- Nextclade data is now checked for compatibility with Nextclade version
- Normalise no longer hardcoded
- Nextclade bump to 1.10.2
- Nextclade data to 2022-02-07T12:00:00Z
- Pangolin bump to 3.1.20

### Fixed
- Pangolin build fix
- NEB VarSkip V2 bed file now corrected
- Nextclade execution options changed as suggested in https://github.com/epi2me-labs/wf-artic/issues/35
- Fixing `--samples` in schema

## [v0.3.10]
### Changes
- `--samples` now changed to `--sample` and `--sample_sheet` introduced
- Schema updated to allow integration with epi2me services
- Added ARTIC V4.1 primer set
- Added NEB VarSkip primer sets
- Updated Nextclade data to `2022-01-05T19:54:31Z`
- Updated Nextclade to `v1.8.0`
- Updated pangolin to `v3.1.17`
- Now outputs all sample VCFs as well as merged VCF
- Pangolin and Nextclade versions can be specified at runtime `--pangolin_version` and `--nextclade_version`
- Pangolin and Nextclade data can be updated at runtime with `--update_data`
- All Nextclade data versions included in repo for user selection with `--nextclade_data_tag`, most recent will be used by default
- Updated max-length of amplicons for NEB long
- Updated NEB-v1a-long bed file to `pool=1` from `pool=0`

## [v0.3.9]
### Changes
- Pangolin update now takes place at build stage, rather than during run

## [v0.3.8]
### Fixed
- Nextclade data was out of date, now calls Omicron

## [v0.3.7]
### Fixed
- Pangolin update step now checks for internet connection

## [v0.3.6]
### Changes
- Pangolin updated to 3.1.16 (Will now call B.1.1.529 successfully)
- Pangolin now auto updates before it runs (internet connectivity required)
- Nextclade updated to 1.5.1

### Fixed
- Segmentation faults in artic-tools vcfcheck code.

## [v0.3.5]
### Added
- Error report if no data at all.
- Telemetry JSON.
- Nextclade errors output as seperate table.
- `--report_detailed` flag to show additional coverage plots
- Parsing of sample sheet to include a `type` column

### Changes
- `sample_name` changed to `sample_id` throughout

## [v0.3.4]
### Added
- Option to add suffix to HTML report name.
- Error message if fastq input file evaluates to null.
- Output Nextflow schema JSON file.
- Output artic JSON file.

## [v0.3.3]
### Changed
- Update nextclade to c++ version 1.3.0, install via bioconda.
- Update aplanat to v0.5.4.

### Added
- V4.1 primer set for spike-seq.
- Tag for pangolin image is now specified in nextflow config.
- Integrate max_softclip_length parameter to be passed into artic minion.
- Output artic.json.

### Fixed
- Parsing of sample_name column from summary files during report curation.

## [v0.3.2]
### Changed
- Updated `fastcat` and `aplanat` versions for standardised software version reporting.

### Fixed
- Empty GVCF file not produced when ARTIC failed.
- `conda` environment file location incorrectly specified in `nextflow.config`

## [v0.3.1]
### Added
- Per-sample bam files now published to output directory.

### Changed
- Data ingress now performed by standard module.

## [v0.3.0]
### Fixed
- Updated medaka to v1.4.3 for model pre-download.
- Work around issue where pyvcf writes QUAL values as '.' and not 0.

### Changed
- Removed the autodetect sample_id option for now.
- Updated default model to be a variant calling one. Although labelled as PromethION specific (`_prom` in name), this model should be preferred on all platforms of non-variant (consensus) platform specific models.
- Derive software versions from CLI rather than conda list.

### Added
- Field `alias` in sample sheet CSV serves as alternative to `sample_name`.
- Added V4 primerscheme to data directory.

## [v0.2.3]
### Changed
- Updated medaka to v1.4.2.
- Updated aplanat to v0.4.0.

## [v0.2.2]
### Added
- Added summary of software parameters section to report.

### Changed
- genotype_variants option can now be used without specifying a path, falling back to the scheme default, if one exists.
- Removed vestigial spike-seq scheme versions.

### Fixed
- Updated allVariants step to normalise REF fields to fix vcf merge issue.
- Prevented nextclade from using all available threads.

## [v0.2.1]
### Fixed
- Intermittent error producing genotyping summary.

## [v0.2.0]
### Added
- Ability to configure depth coverage reporting value.
- Add explicit pins of conda packages.
- Inclusion of SpikeSeq workflow, and reporting.
- Optional auto-detection of sample_id

## [v0.1.4]
### Changed
- Improved display of coverage traces in report.

## [v0.1.3]
### Added
- Check format of sample sheet before executing main workflow.

### Fixed
- Parsing of V1200 .bed file for nextclade report.
- Empty barcode directories are ignored.
- Nextclade report component upgraded to better handle poor data.

## [v0.1.2]
### Fixed
- Recovery after `artic minion` fails.

### Added
- Report item detailing failed analyses.

## [v0.1.1]
### Fixed
- Correct value of wfversion in config.
- Processing of single sample inputs.

## [v0.1.0]
### Added
- Added variant call summary section to report.

## [v0.0.9]
### Changed
- Moved scripting to bin directory.

## [v0.0.8]
### Fixed
- Fix lack of help message when `--help` run.

## [v0.0.7]
### Changed
- Sample sheet is no longer required.
- Sort report items consistently by sample name.
- Nextclade visual will display overlap to primer scheme selected by user.

## [v0.0.6]
### Added
- Support for fragmented amplicons.
- Enabled use of conda profile.

### Changed
- Use custom np-artic package based on 1.3.0-dev branch of original.
- Use nextclade from conda package
- Amended default local executor CPU resource to be more parsimonious.

## [v0.0.5]
### Changed
- Amended report text

## [v0.0.4]
### Changed
- Discretize coverage plots for speed

## [v0.0.3]
### Changed
- Automatically select min/max read lengths base on scheme.

### Added
- Added command-line argument validation.

## [v0.0.2]
* Automation release

### Added
- Continuous deployment.

## [v0.0.1]
* Initial release

### Added
- Basic running of Artic workflow and reporting.

