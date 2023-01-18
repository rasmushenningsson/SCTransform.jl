"""
	sctransform(X::AbstractSparseMatrix, features, params;
                transpose = false,
                feature_id_columns = [:id,:feature_type],
                cell_ind = 1:size(X,2),
                clip=sqrt(size(X,2)/30))

Computes the SCTransform of the sparse matrix `X`, with features as rows and cells as columns.
`features` should be a table (e.g. DataFrame or NamedTuple) with feature annotations.
`params` is a table with reguralized feature parameter estimates - typically estimated using the function `scparams`.

* `transpose` - set to true to transpose the output (cells as rows, features as columns).
* `feature_id_columns` is a vector of column names in `features`. The rows in `features` must be unique based on these columns.
* `cell_ind` is vector of cell indices to include in the output. This is computationally more efficient than subsetting afterwards, but yields the same result.
* `clip` - values less than `-clip` or larger than `clip` are clamped, to reduce the impact of outliers.

See also: [`scparams`](@ref)
"""
function sctransform(X::AbstractSparseMatrix, features, params;
                     transpose = false,
                     feature_id_columns = [:id,:feature_type],
                     cell_ind = 1:size(X,2),
                     clip=sqrt(size(X,2)/30))

	@assert size(X,1)==length(getproperty(features,first(propertynames(features)))) "The number of rows in X and features must match"

	β0 = params.beta0
	β1 = params.beta1
	θ  = params.theta

	# Use features to figure out which rows in params match which rows in X
	# TODO: cleanup code
	param_ids = getproperty(params, feature_id_columns[1])
	feature_ids = getproperty(features, feature_id_columns[1])
	for col in feature_id_columns[2:end]
		param_ids = string.(param_ids, "__<sep>__", getproperty(params, col))
		feature_ids = string.(feature_ids, "__<sep>__", getproperty(features, col))
	end
	feature_ind = indexin(param_ids, feature_ids)
	any(isnothing, feature_ind) && throw(DomainError("Feature ids in `params` does not match ids in `features`."))



	logCellCounts = logcellcounts(X)[cell_ind]

	# TODO: Do not create intermediate X[feature_ind,cell_ind]
	X = X[feature_ind,cell_ind]

	Z = convert(Matrix{Float64}, transpose ? X' : X)

	for j in 1:size(Z,2)
		for i in 1:size(Z,1)
			g,c = transpose ? (j,i) : (i,j)

			μ = exp(β0[g] + β1[g]*logCellCounts[c])
			σ = sqrt(μ + μ^2/θ[g])

			z = (Z[i,j] - μ)/σ
			Z[i,j] = clamp(z, -clip, clip)
		end
	end
	Z
end
