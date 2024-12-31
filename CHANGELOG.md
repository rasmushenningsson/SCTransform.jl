# SCTransform.jl changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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
