function sparsemulsum(sp::SparseVector, d::Vector)::Float64
	@assert length(sp)==length(d)
	s = 0.0
	for (i,v) in zip(SparseArrays.nonzeroinds(sp), nonzeros(sp))
		@inbounds s += v*d[i]
	end
	s
end


"""
	nbnegloglikelihood!(F,G,H,x, y::Vector, c::Vector)

Negative log-likelihood for Negative Binomial distribution, with parameters as in SCTransform.
`y` is the vector of observed counts for the current gene and `c` is the log10 cell counts.
"""
function nbnegloglikelihood!(F,G,H,x, y::Vector, c::Vector)
	r,β0,β1 = x[1],x[2],x[3]
	rSign = sign(r)
	r = abs(r)
	n = length(y)

	logμ = β0.+β1.*c
	μ = exp.(logμ)
	logrμ = log.(r.+μ) # needed by F and G
	

	if G != nothing
		# compute gradient and store in G
		G[1] = rSign*-(n*(-digamma(r) + log(r) + 1) + sum( digamma.(r.+y) .- logrμ .- (r.+y)./(r.+μ)))
		G[2] = -(sum( -(r.+y).*μ./(r.+μ) .+ y ))
		G[3] = -(sum( -(r.+y).*c.*μ./(r.+μ) .+ y.*c ))
	end
	if H != nothing
		# compute Hessian and store in H
		dldr2 = -(n*(-trigamma(r) + 1/r) + sum(trigamma.(r.+y) .- 2.0./(r.+μ) .+ (r.+y)./(r.+μ).^2) )
		dldβ02 = sum((r.+y).*( μ./(r.+μ) .- μ.^2 ./ (r.+μ).^2))
		dldβ12 = sum(c.*(r.+y) .* (μ.*c./(r.+μ) .- c.*μ.^2 ./(r.+μ).^2))
		dldrdβ0 = sum(μ .* ( 1.0./(r.+μ) .- (r.+y)./(r.+μ).^2 ))
		dldrdβ1 = sum(μ.*c .* ( 1.0./(r.+μ) .- (r.+y)./(r.+μ).^2 ))
		dldβ0dβ1 = sum((r.+y).*c .* ( μ./(r.+μ) .- μ.^2 ./(r.+μ).^2 ))

		H[1,1] = dldr2
		H[2,2] = dldβ02
		H[3,3] = dldβ12
		H[1,2] = H[2,1] = rSign*dldrdβ0
		H[1,3] = H[3,1] = rSign*dldrdβ1
		H[2,3] = H[3,2] = dldβ0dβ1
	end
	if F != nothing
		return -(n*(-loggamma(r) + r*log(r)) + sum( loggamma.(r.+y) .- r.*logrμ .+ y.*logμ .- y.*logrμ ))
	end
	nothing
end


"""
	poissonnegloglikelihood!(F,G,H,x, y::AbstractVector, c::Vector)

Negative log-likelihood for Poisson distribution, with parameters as in SCTransform.
y is the vector of observed counts for the current gene and `c` is the log10 cell counts.
"""
function poissonnegloglikelihood! end

function poissonnegloglikelihood!(F,G,H,x, y::AbstractVector, c::Vector)
	β0,β1 = x[1],x[2]

	logμ = β0.+β1.*c
	μ = exp.(logμ)

	if G != nothing
		# compute gradient and store in G
		G[1] = -sum(y .- μ)
		G[2] = -sum(c.*(y.-μ))
	end
	if H != nothing
		# compute Hessian and store in H
		dldβ02   = sum(μ)
		dldβ12   = sum(c.^2.0.*μ)
		dldβ0dβ1 = sum(c.*μ)

		H[1,1] = dldβ02
		H[2,2] = dldβ12
		H[1,2] = H[2,1] = dldβ0dβ1
	end
	if F != nothing
		return -sum(y.*logμ .- μ)
	end
	nothing
end

function poissonnegloglikelihood!(F,G,H,x, y::SparseVector, c::Vector)
	β0,β1 = x[1],x[2]
	N = length(y)

	nzInd = SparseArrays.nonzeroinds(y)
	nzVal = nonzeros(y)

	sF = 0.0
	sG1 = sG2 = 0.0
	sH11 = sH12 = sH22 = 0.0

	sparseInd = 1
	nextNZInd = isempty(nzInd) ? 0 : nzInd[1]
	for i=1:N
		@inbounds ci = c[i]
		μi = exp(β0 + β1*ci)
		H != nothing && (sH11 += μi; sH12 += ci*μi; sH22 += ci*ci*μi)
		if i!=nextNZInd # yᵢ=0
			F != nothing && (sF += μi)
			G != nothing && (sG1 += μi; sG2 += ci*μi)
		else # yᵢ≠0
			@inbounds yi = nzVal[sparseInd]
			F != nothing && (sF += μi - β0*yi - β1*ci*yi)
			G != nothing && (sG1 += μi-yi; sG2 += ci*(μi-yi))
			sparseInd += 1
			sparseInd<=length(nzInd) && (@inbounds nextNZInd = nzInd[sparseInd])
		end
	end
	G != nothing && (G[1]=sG1; G[2]=sG2)
	H != nothing && (H[1,1]=sH11; H[1,2]=H[2,1]=sH12; H[2,2]=sH22)
	F != nothing && return sF
	nothing
end
# fast but less accurate
# function poissonnegloglikelihood!(F,G,H,x, y::SparseVector, c::Vector)
# 	β0,β1 = x[1],x[2]

# 	μSum = 0.0
# 	cμSum = 0.0
# 	c2μSum = 0.0
# 	for ci in c
# 		logμi = β0 + β1*ci
# 		μi = exp(logμi)
# 		μSum += μi
# 		cμSum += ci*μi
# 		c2μSum += ci*ci*μi
# 	end
# 	ySum  = sum(y)
# 	ycSum = sparsemulsum(y,c)
# 	if G != nothing
# 		G[1] = μSum - ySum
# 		G[2] = cμSum - ycSum
# 	end
# 	if H != nothing
# 		H[1,1] = μSum
# 		H[2,2] = c2μSum
# 		H[1,2] = H[2,1] = cμSum
# 	end
# 	if F != nothing
# 		return μSum - β0*ySum - β1*ycSum # sum(y.*(β0.+β1.*c)) = β0*sum(y) + β1*sum(y.*c)
# 	end
# 	nothing
# end


"""
	thetanegloglikelihood!(F,G,H,x, y::AbstractVector, μ::Vector, logμ::Vector)

Negative log-likelihood for Negative Binomial distribution, with fixed `μ` parameter. Useful for estimating `θ` after `μ` has been estimated in a Poisson model.
`y` is the vector of observed counts for the current gene.
"""
function thetanegloglikelihood! end

function thetanegloglikelihood!(F,G,H,x, y::AbstractVector, μ::Vector, logμ::Vector)
	r = x[1]
	rSign = sign(r)
	r = abs(r)
	n = length(y)

	logrμ = log.(r.+μ) # needed by F and G
	
	if G != nothing
		# compute gradient and store in G
		G[1] = rSign*-(n*(-digamma(r) + log(r) + 1) + sum( digamma.(r.+y) .- logrμ .- (r.+y)./(r.+μ)))
	end
	if H != nothing
		# compute Hessian and store in H
		dldr2 = -(n*(-trigamma(r) + 1/r) + sum(trigamma.(r.+y) .- 2.0./(r.+μ) .+ (r.+y)./(r.+μ).^2) )
		H[1,1] = dldr2
	end
	if F != nothing
		return -(n*(-loggamma(r) + r*log(r)) + sum( loggamma.(r.+y) .- r.*logrμ .+ y.*logμ .- y.*logrμ ))
	end
	nothing
end

function thetanegloglikelihood!(F,G,H,x, y::SparseVector, μ::Vector, logμ::Vector)
	@assert length(y)==length(μ)==length(logμ)
	r = x[1]
	rSign = sign(r)
	r = abs(r)
	N = length(y)

	nzInd = SparseArrays.nonzeroinds(y)
	nzVal = nonzeros(y)

	sF = 0.0
	sG = 0.0
	sH = 0.0

	loggammar = F != nothing ? loggamma(r) : 0.0
	digammar  = G != nothing ? digamma(r)  : 0.0
	trigammar = H != nothing ? trigamma(r) : 0.0

	logr = log(r)

	sparseInd = 1
	nextNZInd = isempty(nzInd) ? 0 : nzInd[1]
	for i=1:N
		if i!=nextNZInd # yᵢ=0
			@inbounds μi = μ[i]
			local logrμi::Float64
			(F != nothing || G != nothing) && (logrμi = log(r+μi))
			F != nothing && (sF += r*logr - r*logrμi)
			rμRep = 1.0/(r+μi)
			G != nothing && (sG += 1 + logr - logrμi - r*rμRep)
			H != nothing && (sH += 1.0/r - 2.0*rμRep + r*rμRep*rμRep)
		else # yᵢ≠0
			@inbounds μi = μ[i]
			local logμi::Float64
			F != nothing && (@inbounds logμi = logμ[i])
			(F != nothing || G != nothing) && (logrμi = log(r+μi))

			@inbounds yi = nzVal[sparseInd]
			rμRep = 1.0/(r+μi)
			F != nothing && (sF += loggamma(r+yi) - loggammar + r*logr + yi*logμi - (r+yi)*logrμi)
			G != nothing && (sG += digamma(r+yi) - digammar + 1 + logr - logrμi - (r+yi)*rμRep)
			H != nothing && (sH += trigamma(r+yi) - trigammar + 1.0/r - 2.0*rμRep + (r+yi)*rμRep*rμRep)
			sparseInd += 1
			# nextNZInd = sparseInd>length(nzInd) ? 0 : nzInd[sparseInd]
			sparseInd<=length(nzInd) && (@inbounds nextNZInd = nzInd[sparseInd])
		end
	end
	G != nothing && (G[1] = -rSign*sG)
	H != nothing && (H[1,1] = -sH)
	F != nothing && return -sF
	nothing
end
# fast but less accurate
# function thetanegloglikelihood!(F,G,H,x, y::SparseVector, μ::Vector, logμ::Vector)
# 	@assert length(y)==length(μ)==length(logμ)
# 	r = x[1]
# 	rSign = sign(r)
# 	r = abs(r)
# 	n = length(y)

# 	logrμSum = 0.0
# 	rμRepSum = 0.0
# 	rμ2RepSum = 0.0
# 	for μi in μ
# 		logrμSum += log(r+μi)
# 		rμRep = 1.0/(r+μi)
# 		rμRepSum += rμRep
# 		rμ2RepSum += rμRep*rμRep
# 	end

# 	if G != nothing
# 		s = (n-nnz(y))*digamma(r) - r*rμRepSum
# 		for (i,yi) in zip(SparseArrays.nonzeroinds(y), nonzeros(y))
# 			@inbounds s += digamma(r+yi) - yi/(r+μ[i])
# 		end
# 		G[1] = rSign*-(n*(-digamma(r) + log(r) + 1) + s - logrμSum)
# 	end
# 	if H != nothing
# 		s = (n-nnz(y))*trigamma(r) - 2.0*rμRepSum + r*rμ2RepSum
# 		for (i,yi) in zip(SparseArrays.nonzeroinds(y), nonzeros(y))
# 			@inbounds s += trigamma(r+yi) + yi/(r+μ[i])^2
# 		end
# 		H[1,1] = -(n*(-trigamma(r) + 1/r) + s)
# 	end
# 	if F != nothing
# 		s = (n-nnz(y))*loggamma(r)
# 		for (i,yi) in zip(SparseArrays.nonzeroinds(y), nonzeros(y))
# 			@inbounds s += loggamma(r+yi) + yi*logμ[i] - yi*log(r+μ[i])
# 		end
# 		return -(n*(-loggamma(r) + r*log(r)) + s - r*logrμSum)
# 	end
# 	nothing
# end



"""
	nbparams(y::Vector, c::Vector; verbose=false)

Fits a Negative Binomial distribution to the given data, with parameters `μ = exp(β₀+β₁c)`
and dispersion parameter `θ` by performing ml-estimatation of `θ`,`β₀`,`β₁`.
`y` is the vector of observed counts for the current gene and `c` is the log10 cell counts.

See `nbparamsbypoisson` for a faster, but less accurate, alternative.
"""
function nbparams(y::Vector, c::Vector; verbose=false)
	res = Optim.optimize(Optim.only_fgh!((F,G,H,x)->nbnegloglikelihood!(F,G,H,x,y,c)), [0.1, 0., 0.], Optim.Newton())
	verbose && @info res
	verbose && !Optim.converged(res) && @warn "Negative Binomial ml-estimation of θ,β0,β1 failed."
	θ,β0,β1 = Optim.minimizer(res)
	abs(θ),β0,β1
end

struct NBByPoissionScratch
	βinit::Vector{Float64}
	θinit::Vector{Float64}
	μ::Vector{Float64}
	logμ::Vector{Float64}
end
NBByPoissionScratch(N::Int) = NBByPoissionScratch(zeros(2),zeros(1),zeros(N),zeros(N))

"""
	nbparamsbypoisson(y::Union{Vector,SparseVector}, c::Vector; verbose=false)

First fits a Poisson model to the given data, with parameter `μ = exp(β₀+β₁c)` by ml-estimation of `β₀`,`β₁`.
Then computes the ml-estimate `θ` of a Negative Binomial model with `μ` fixed.
`y` is the vector of observed counts for the current gene and `c` is the log10 cell counts.

See `nbparams` for a slower, but more accurate, alternative.
"""
function nbparamsbypoisson(y::AbstractVector, c::Vector; verbose=false, scratch=NBByPoissionScratch(length(y)))
	scratch.βinit .= 0.0
	res = Optim.optimize(Optim.only_fgh!((F,G,H,x)->poissonnegloglikelihood!(F,G,H,x,y,c)), scratch.βinit, Optim.Newton())

	verbose && @info res
	verbose && !Optim.converged(res) && @warn "Poisson ml-estimation of β0,β1 failed."
	β0,β1 = Optim.minimizer(res)

	# estimate θ given β₀,β₁
	scratch.logμ .= β0.+β1.*c
	scratch.μ    .= exp.(scratch.logμ)
	scratch.θinit[1] = 0.1
	res2 = Optim.optimize(Optim.only_fgh!((F,G,H,x)->thetanegloglikelihood!(F,G,H,x,y,scratch.μ,scratch.logμ)), scratch.θinit, Optim.Newton())

	verbose && @info res2
	verbose && !Optim.converged(res2) && @warn "ML-estimation of θ failed."
	θ = abs(Optim.minimizer(res2)[1])

	θ,β0,β1
end



function thetastandarderror(y::AbstractVector,c::Vector,θ::Real,β0::Real,β1::Real; scratch=NBByPoissionScratch(length(y)))
	scratch.logμ .= β0.+β1.*c
	scratch.μ    .= exp.(scratch.logμ)
	H = zeros(1,1)
	thetanegloglikelihood!(nothing,nothing,H,θ,y,scratch.μ,scratch.logμ)
	1/sqrt(max(0.0,H[1,1]))
end

function thetaconfidenceinterval(y::AbstractVector,c::Vector,θ::Real,β0::Real,β1::Real; kwargs...)
	pm = 1.9599639845400576 * thetastandarderror(y,c,θ,β0,β1; kwargs...)
	(θ-pm,θ+pm)
end
