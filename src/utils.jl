"""
	splitrange(r::UnitRange, nparts::Integer)

Splits a range in `nparts` number of parts of equal length.
"""
function splitrange(r::UnitRange{T}, nbrParts::Integer) where T<:Real
	s = first(r)
	d,r = divrem(length(r),nbrParts)
	out = Vector{UnitRange{T}}(undef, nbrParts)
	for i=1:nbrParts
		len = d+(i<=r)
		out[i] = range(s, length=len)
		s += len
	end
	out
end
