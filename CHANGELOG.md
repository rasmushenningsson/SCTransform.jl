# SCTransform.jl changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Breaking

WIP
* `scparams` now uses `feature_mask` when computing `logcellcounts`. This is important there are multiple modalities in the data (e.g. Gene Expression and Antibody counts.)
* Caching of `scparams` results now uses `StableHashTraits` for hashing.
* `sctransform` now takes `feature_mask` parameter, which controls which features are used to compute `logcellcounts`.