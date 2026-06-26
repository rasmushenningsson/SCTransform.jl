# SCTransform.jl changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.4.1] - 2026-06-25

### Fixed

* Loosen compat with Optim (to 1.2,2), and NLSolversBase to (7,8).

## [0.4.0] - 2026-06-25

### Breaking

* Features for which the parameter estimation fails (very rare) are now marked as outliers instead of being removed. For this reason, previously cached parameters will be recomputed.

### Added

* Low-level driver function `compute_scparams`.

### Fixed

* Fix potential (but unlikely) data race due to a BitVector being accessed from multiple threads.
* Failed inference is now report with one warning instead of one per failed feature.

### Changed

* Internal refactoring to make it work better with the upcoming SingleCellProjections.jl v0.5.
* Require at least Julia 1.10 (current LTS).

## [0.3.1] - 2024-12-31

### Fixed

* Compat with StableHashTraits v2.

## [0.3.0] - 2024-11-08

### Breaking

* `stable_hash` will now use hash version 4 by default, since that is the only one compatible with Julia 1.11+. It will not affect any results, but any sctransform parameters that were cached on disk will be recomputed.

## [0.2.0] - 2024-06-18

### Breaking

* `scparams` now uses `feature_mask` when computing `logcellcounts`. This is important if there are multiple modalities in the data (e.g. Gene Expression and Antibody counts.)
* `sctransform` now takes `feature_mask` parameter, which controls which features are used to compute `logcellcounts`. Defaults to only using "Gene Expression" features.
* Caching of `scparams` results now uses `StableHashTraits` for hashing.
