# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
- Updated `fastcat` and `aplanat` versions for standardised software version
  reporting.
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
- Updated default model to be a variant calling one. Although labelled as
  PromethION specific (`_prom` in name), this model should be preferred
  on all platforms of non-variant (consensus) platform specific models.
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
- genotype_variants option can now be used without specifying a path, falling
  back to the scheme default, if one exists.
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
###
- Added command-line argument validation.

## [v0.0.2]

Automation release

### Added
- Continuous deployment.


## [v0.0.1]

Initial release

### Added
- Basic running of Artic workflow and reporting.
