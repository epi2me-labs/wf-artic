# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
