module SCTransform

export
	scparams,
	scparams_estimate,
	scparams_detect_outliers,
	scparams_bandwidth,
	scparams_regularize,
	sctransform

using LinearAlgebra
using SparseArrays
using Statistics
using SpecialFunctions
using Optim
using KernelDensitySJ
using ProgressMeter
using Scratch
using SHA
using DelimitedFiles
using CodecZlib

include("utils.jl")
include("table.jl")
include("nbregression.jl")
include("smooth.jl")
include("cache.jl")
include("params.jl")
include("transform.jl")


const scparams_scratch_space = Ref{String}()

function __init__()
    global scparams_scratch_space[] = @get_scratch!("SCParams")
end

end
