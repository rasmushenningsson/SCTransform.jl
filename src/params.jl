

# Assumes each column is a gene
# logcellcounts(X::SparseMatrixCSC) = log10.(max.(1,vec(sum(X;dims=2))))

# Assumes each column is a cell
logcellcounts(X::SparseMatrixCSC) = log10.(max.(1,vec(sum(X;dims=1))))

# Assumes each column is a gene
# function loggenemean(X::SparseMatrixCSC)
# 	N,P = size(X)
# 	lgm = zeros(P)
# 	V = nonzeros(X)
# 	for j=1:P
# 		s = 0.0
# 		for k in nzrange(X,j)
# 			s += log(V[k]+1)
# 		end
# 		lgm[j] = s/N
# 	end
# 	log10.(exp.(lgm).-1)
# end


# Assumes each column is a cell
function loggenemean(X::SparseMatrixCSC)
	N = size(X,2)
	lgm = sum(log1p, X; dims=2) / N
	log10.(exp.(lgm).-1)
end



function scparams_poisson_worker(channel, progress, N, feature_mask, feature_names, logCellCounts, min_cells, θ, β0, β1, θSE, kept)
	scratch = NBByPoissionScratch(N)

	while true
		item = take!(channel)
		isnothing(item) && break # no more chunks to process


		X,feature_offset = item
		for j in 1:size(X,2)
			if length(nzrange(X,j)) >= min_cells && feature_mask[j+feature_offset]
				sparseY = X[:,j]

				j2 = j+feature_offset
				try
					θ[j2],β0[j2],β1[j2] = nbparamsbypoisson(sparseY,logCellCounts,scratch=scratch)
					θSE[j2] = thetastandarderror(sparseY,logCellCounts,θ[j2],β0[j2],β1[j2],scratch=scratch)
					kept[j2] = true
				catch e
					@warn "Failed to compute SCTransform parameters for feature $(feature_names[j2]), skipping."
				end
			end
			isnothing(progress) || next!(progress)
		end
		yield() # is this needed when we are using a channel?
	end
end


function scparams_nb_worker(channel, progress, N, feature_mask, feature_names, logCellCounts, min_cells, θ, β0, β1, θSE, kept)
	while true
		item = take!(channel)
		isnothing(item) && break # no more chunks to process

		X,feature_offset = item
		for j in 1:size(X,2)
			if length(nzrange(X,j)) >= min_cells && feature_mask[j+feature_offset]
				sparseY = X[:,j]

				j2 = j+feature_offset
				try
					y = convert(Vector,sparseY) # TODO: get rid of conversion to full vector by fixing sparse version of nbparams
					θ[j2],β0[j2],β1[j2] = nbparams(y,logCellCounts)
					θSE[j2] = thetastandarderror(sparseY,logCellCounts,θ[j2],β0[j2],β1[j2])
					kept[j2] = true
				catch e
					@warn "Failed to compute SCTransform parameters for feature $(feature_names[j2]), skipping."
				end
			end
			isnothing(progress) || next!(progress)
		end
		yield() # is this needed when we are using a channel?
	end
end



function scparams_estimate(::Type{T}, X::AbstractSparseMatrix{Tv,Ti};
                           method=:poisson,
                           min_cells::Integer=5,
                           feature_mask = trues(size(X,1)),
                           feature_names,
                           chunk_size=100,
                           nthreads=Threads.nthreads(),
                           channel_size=nthreads*4,
                           verbose=true,
                          ) where {T,Tv<:Real,Ti<:Integer}
	nthreads = max(nthreads,1)
	P,N = size(X)

	logCellCounts=logcellcounts(X)
	logGeneMean=loggenemean(X)

	@assert method in (:poisson, :nb) "Method must be :poisson or :nb"

	@assert length(feature_mask) == P

	θ    = zeros(P)
	β0   = zeros(P)
	β1   = zeros(P)
	θSE  = zeros(P)
	kept = falses(P)

	progress = verbose ? Progress(P; desc="Estimating sc parameters: ") : nothing


	channel = Channel{Union{Nothing,Tuple{SparseMatrixCSC{Tv,Ti},Int}}}(channel_size)


	workers = map(1:nthreads) do _
		if method==:poisson
			Threads.@spawn scparams_poisson_worker(channel, progress,
			                                       N, feature_mask, feature_names,
			                                       logCellCounts, min_cells,
			                                       θ, β0, β1, θSE, kept)
		elseif method==:nb
			Threads.@spawn scparams_nb_worker(channel, progress,
			                                  N, feature_mask, feature_names,
			                                  logCellCounts, min_cells,
			                                  θ, β0, β1, θSE, kept)
		end
	end



	colptr_curr = first.(nzrange.(Ref(X),1:N))
	colptr_end = last.(nzrange.(Ref(X),1:N))
	rowval = rowvals(X)
	nzval = nonzeros(X)

	rowval_scratch = Vector{Ti}() # will grow but get reused between chunks
	nzval_scratch  = Vector{Tv}() # will grow but get reused between chunks

	for feature_range in Iterators.partition(1:P, chunk_size)
		if any(feature_mask[feature_range])
			colptr_chunk = Vector{Ti}(undef, N+1)

			for j in 1:N
				colptr_chunk[j] = length(rowval_scratch)+1

				c = colptr_curr[j]

				while c<=colptr_end[j] && rowval[c]<=last(feature_range)
					push!(rowval_scratch, rowval[c] - first(feature_range) + 1)
					push!(nzval_scratch, nzval[c])
					c += 1
				end

				colptr_curr[j] = c
			end
			colptr_chunk[end] = length(rowval_scratch)+1

			rowval_chunk = copy(rowval_scratch)
			nzval_chunk = copy(nzval_scratch)

			empty!(rowval_scratch)
			empty!(nzval_scratch)

			chunk = SparseMatrixCSC(length(feature_range), N, colptr_chunk, rowval_chunk, nzval_chunk)
			chunk = permutedims(chunk,(2,1)) # transpose

			push!(channel, (chunk,first(feature_range)-1))
		else # no features used in this chunk, just step colptr forward
			# TODO: Do a binary search instead?
			#       Or better, skip multiple chunks at the same time if possible. (Lazy skipping + binary search when we need to.)
			for j in 1:N
				c = colptr_curr[j]
				while c<=colptr_end[j] && rowval[c]<=last(feature_range)
					c += 1
				end
				colptr_curr[j] = c
			end
		end
	end

	# Tell workers to stop
	for i in 1:nthreads
		push!(channel, nothing)
	end

	wait.(workers)
	isnothing(progress) || finish!(progress)

	@assert any(kept) "SC Parameter estimation failed - no features remaining."

	T((; featureInd=(1:P)[kept], logGeneMean=logGeneMean[kept], beta0_estimate=β0[kept], beta1_estimate=β1[kept], theta_estimate=θ[kept], thetaSE_estimate=θSE[kept]))
end




"""
windowmedian(x::AbstractVector, y::AbstractVector, d::Real)

Simple implementation of windowed median.
`x` must be sorted.
For element `i` in `y`, the median is computed for all indices with x coordinate in the interval `x[i]±d/2`.
"""
function windowedmedian(x::AbstractVector, y::AbstractVector, d::Real)
	# TODO: Optimize using heaps
	@assert issorted(x)
	@assert length(x)==length(y)
	@assert d>=0

	N = length(x)
	m = zeros(N)
	iStart = 1
	iEnd = 2
	for i=1:length(x)
		while x[i]-x[iStart] > d/2
			iStart += 1
		end
		while iEnd<=N && x[iEnd]-x[i] <= d/2
			iEnd += 1
		end
		m[i] = median(view(y,iStart:iEnd-1))
	end
	m
end

function windowedoutlierscore(x::AbstractVector, y::AbstractVector, d::Real)
	med = windowedmedian(x, y, d)
	dev = y.-med
	mad = windowedmedian(x, abs.(dev), d).*1.4826022185056018 # commonly used scale factor 1/Φ⁻¹(3/4)
	abs.(dev./(mad))
end


function scparams_detect_outliers(::Type{T}, params) where T
	mask = params.thetaSE_estimate.<=4 # flag extremely unreliable theta estimates as outliers
	rows = (1:length(params.logGeneMean))[mask] # to match original rows later

	# Work with sorted data
	perm = sortperm(params.logGeneMean[mask])
	logGeneMean  = permute!(params.logGeneMean[mask], perm)
	beta0_estimate = permute!(params.beta0_estimate[mask], perm)
	beta1_estimate = permute!(params.beta1_estimate[mask], perm)
	theta_estimate = permute!(params.theta_estimate[mask], perm)

	windowSize = (logGeneMean[end]-logGeneMean[1]) / sqrt(length(logGeneMean))

	outliers_beta0 = windowedoutlierscore(logGeneMean, beta0_estimate, windowSize) .> 10
	outliers_beta1 = windowedoutlierscore(logGeneMean, beta1_estimate, windowSize) .> 10
	outliers_theta = windowedoutlierscore(logGeneMean, theta_estimate, windowSize) .> 10
	outliers = outliers_beta0 .| outliers_beta1 .| outliers_theta

	# Go back to original order
	invpermute!(outliers, perm)
	outliers2 = trues(length(params.featureInd))
	outliers2[rows] .= outliers

	T(hcat_tables(params, (;outlier=outliers2)))
end

scparams_detect_outliers(params::T) where T = scparams_detect_outliers(T,params)
scparams_detect_outliers(params::NamedTuple) = scparams_detect_outliers(NamedTuple,params) # strip NamedTuple parameters



function scparams_bandwidth(params; kwargs...)
	bwsj(params.logGeneMean[.!params.outlier]; kwargs...)
end


function scparams_regularize(::Type{T}, params, bw) where T
	mask = .!params.outlier

	logGeneMean  = params.logGeneMean[mask]
	beta0_estimate = params.beta0_estimate[mask]
	beta1_estimate = params.beta1_estimate[mask]
	theta_estimate = params.theta_estimate[mask]


	smoothed_beta0 =       gaussiansmoothing(logGeneMean,        beta0_estimate,  bw, params.logGeneMean)
	smoothed_beta1 =       gaussiansmoothing(logGeneMean,        beta1_estimate,  bw, params.logGeneMean)
	smoothed_theta = 10.0.^gaussiansmoothing(logGeneMean, log10.(theta_estimate), bw, params.logGeneMean)

	t1 = remove_columns(params, (:beta0_estimate, :beta1_estimate, :theta_estimate, :thetaSE_estimate))
	t2 = (; beta0=smoothed_beta0,
	        beta1=smoothed_beta1,
	        theta=smoothed_theta)
	T(hcat_tables(t1, t2))
end

scparams_regularize(params::T, bw) where T = scparams_regularize(T,params,bw)
scparams_regularize(params::NamedTuple, bw) = scparams_regularize(NamedTuple,params,bw) # strip NamedTuple parameters




function scparams(::Type{T}, X::AbstractSparseMatrix, features;
                  method=:poisson,
                  min_cells::Integer=5,
                  feature_type = hasproperty(features, :feature_type) ? "Gene Expression" : nothing,
                  feature_mask = feature_type !== nothing ? features.feature_type.==feature_type : trues(size(X,1)),
                  feature_names = hasproperty(features,:name) ? features.name : features.id,
                  use_cache=true,
                  cache_read=use_cache,
                  cache_write=use_cache,
                  verbose=true,
                  kwargs...) where T
	P,N = size(X)
	length(feature_names) == P || throw(DimensionMismatch("The number of rows in the count matrix and the number of features do not match."))

	if cache_read || cache_write
		h = _scparams_checksum(X,method,min_cells,feature_mask)
		fn_cached = joinpath(scparams_scratch_space[],string(h,".tsv.gz"))
	end

	# check if the results is already in the cache
	params = nothing
	if cache_read && isfile(fn_cached)
		params = _scparams_cache_load(fn_cached, P, N, method, min_cells)
		if params !== nothing
			@info "SCTransform parameters loaded from cache."
			touch(fn_cached) # update file timestamp (in case we want to remove old cached files later)
		end
	end

	if params === nothing
		# 1. fit per gene
		params = scparams_estimate(NamedTuple, X; method, min_cells, feature_mask, feature_names, verbose, kwargs...)

		# 2. detect outliers
		params = scparams_detect_outliers(NamedTuple, params)

		# 3. Compute bandwidth, ignoring outlier values
		# NB: In the SCTransform paper, this was multiplied by 3. However, because of different scaling in ksmooth, a factor of quantile(Normal(),0.75)/0.25 ≈ 2.6979590007843273 would have been correct. Since gaussiansmoothing below uses the same scaling as bwsj, we do not rescale.
		bw = scparams_bandwidth(params) # TODO: pass on bwsj kwargs?
		verbose && @info "Bandwidth $bw"

		# 4. Evaluate smoothed function at all genes
		params = scparams_regularize(NamedTuple, params, bw)

		if cache_write
			_scparams_cache_save(fn_cached, params, P, N, method, min_cells)
		end
	end

	# 5. Combine with feature annotations
	features_subset = subset_rows(features,params.featureInd)
	params = remove_columns(params, (:featureInd,))
	T(hcat_tables(features_subset, params))
end

function scparams(X::AbstractSparseMatrix, features::T; kwargs...) where T
	scparams(T, X, features; kwargs...)
end
function scparams(X::AbstractSparseMatrix, features::NamedTuple; kwargs...)
	scparams(NamedTuple, X, features; kwargs...) # strip NamedTuple parameters
end

