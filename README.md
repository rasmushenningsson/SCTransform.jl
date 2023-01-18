# SCTransform.jl

This is a Julia implementation of [sctransform](https://github.com/satijalab/sctransform), based on the paper [Normalization and variance stabilization of single-cell RNA-seq data using regularized negative binomial regression](https://doi.org/10.1186/s13059-019-1874-1) by Hafemaister and Satija.

This is an independent implementation, the original authors were not involved.
Here is [the original R implementation of sctransform](https://github.com/satijalab/sctransform).

## Usage
The doc strings for `scparams` and `sctransform` describe the usage.
Also see the [SingleCellProjections.jl](https://github.com/rasmushenningsson/SingleCellProjections.jl) package which provides a framework for working with Single Cell expression data and supports `SCTransform.jl`.

Note that [threading](https://docs.julialang.org/en/v1/manual/multi-threading/) should be enabled in Julia for `SCTransform.jl` to work optimally.
