const SCPARAMS_VERSION = v"0.2.0"

function _scparams_checksum(X::SparseMatrixCSC, method::Symbol, min_cells::Int, feature_mask::BitVector)
	@assert method in (:poisson, :nb) "Method must be :poisson or :nb"
	P,N = size(X)

	tup = (P, N, X.colptr, X.rowval, X.nzval, method, min_cells, feature_mask, SCPARAMS_VERSION)
	bytes2hex(stable_hash(tup; version=4))
end

function _scparams_cache_convert_column(data, type)
	if type=="Float64"
		parse.(Float64, data)
	elseif type=="Int64"
		parse.(Int64, data)
	elseif type=="Bool"
		parse.(Bool, data)
	elseif type=="Float32"
		parse.(Float32, data)
	elseif type=="Int32"
		parse.(Int32, data)
	else
		error("Unknown column type: $t")
	end
end

function _scparams_cache_load(fn, P, N, method, min_cells)
	try
		open(fn, "r") do io_raw
			io = GzipDecompressorStream(io_raw)

			line = readline(io)
			@assert startswith(line, "#version=")
			version = parse(VersionNumber, line[10:end])
			@assert version==SCPARAMS_VERSION

			line = readline(io)
			@assert startswith(line, "#P=")
			cached_P = parse(Int, line[4:end])
			@assert P==cached_P

			line = readline(io)
			@assert startswith(line, "#N=")
			cached_N = parse(Int, line[4:end])
			@assert N==cached_N

			line = readline(io)
			@assert startswith(line, "#method=")
			cached_method = Symbol(line[9:end])
			@assert method==cached_method

			line = readline(io)
			@assert startswith(line, "#min_cells=")
			cached_min_cells = parse(Int, line[12:end])
			@assert min_cells==cached_min_cells

			line = readline(io)
			@assert startswith(line, "#")
			types = split(line[2:end],'\t')

			data, header = readdlm(io, '\t', String; header=true)
			header = Symbol.(header)

			close(io)

			@assert length(types)==length(header)
			columns = [_scparams_cache_convert_column(data[:,i],t) for (i,t) in enumerate(types)]

			NamedTuple(h=>c for (h,c) in zip(header,columns))
		end
	catch ex
		@warn "Failed to load cached SCTransform parameters from cache."
		showerror(stdout, ex, catch_backtrace())
		nothing
	end
end

function _scparams_cache_save(fn, params, P, N, method, min_cells)
	open(fn, "w") do io_raw
		io = GzipCompressorStream(io_raw)
		println(io, "#version=", SCPARAMS_VERSION)
		println(io, "#P=", P)
		println(io, "#N=", N)
		println(io, "#method=", method)
		println(io, "#min_cells=", min_cells)

		types = [eltype(getproperty(params,col)) for col in propertynames(params)]
		print(io, '#')
		join(io, types, '\t')
		println(io)

		join(io, propertynames(params), '\t')
		println(io)

		# X = string.(Matrix(params))
		data = [string.(getproperty(params,col)) for col in propertynames(params)]
		data = reduce(hcat,data)

		writedlm(io, data, '\t')
		close(io)
	end
end
