# simple implementation, should be optimized to avoid n² complexity!
function gaussiansmoothing(x::AbstractVector{<:T1}, y::AbstractVector{<:S}, bandwidth::T2, xeval::AbstractVector{<:T3}) where {T1<:Real,S<:Real,T2<:Real,T3<:Real}
	n = length(x)
	m = length(xeval)
	@assert length(y)==n
	c = 1/(2*bandwidth^2)
	smoothed = similar(xeval, promote_type(S,T1,T2,T3))
	w = similar(x, promote_type(T1,T2,T3))
	for i=1:m
		w .= exp.(-c*(x.-xeval[i]).^2)
		smoothed[i] = y'w / sum(w)
	end
	smoothed
end
