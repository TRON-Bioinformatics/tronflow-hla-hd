# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

### Changed

### Removed


## [0.3.0] - 2024-03-25

### Changed

- Replaced GATKs SamToFastq by bedtools bamtofastq to avoid issue with unmated reads


## [0.2.0] - 2024-03-25

### Added

- Add support for BAM files as input files. It extracts the reads aligned to the MHC region in the canonical 
chromosome 6 + all reads from chromosome 6 alt contigs + unmapped reads and converts into FASTQs


## [0.1.0] - 2022-11-02

- First release!


[unreleased]: https://github.com/tron-bioinformatics/tronflow-hla-hd/compare/v0.3.0...HEAD
[0.3.0]: https://github.com/tron-bioinformatics/tronflow-hla-hd/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/tron-bioinformatics/tronflow-hla-hd/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/tron-bioinformatics/tronflow-hla-hd/releases/tag/v0.1.0