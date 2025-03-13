module SCTransform

export
	scparams,
	sctransform


# Use public keyword in Julia versions where it is available
if VERSION >= v"1.11.0-DEV.469"
    let str = """
        public scparams_impl
        """
        eval(Meta.parse(str))
    end
end

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
using StableHashTraits

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
