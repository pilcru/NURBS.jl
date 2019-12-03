module Bases

import IterTools: groupby, imap
import LinearAlgebra: dot
import ..Utils: Interval

export Basis, Basis1D, domain, deriv, order,
       BSplineBasis,
       NURBSBasis


# Abstract Basis types
# ========================================================================

abstract type Basis end
abstract type Basis1D <: Basis end

mutable struct BasisFunction1D{B<:Basis1D}
    basis::B
    index::Int
    deriv::Int

    function BasisFunction1D{B}(basis, index, deriv) where {B<:Basis1D}
        if !(1 <= index <= length(basis)) throw(BoundsError()) end
        if !(0 <= deriv <= nderivs(basis)) throw(ArgumentError("Differentiation order not supported")) end
        new(basis, index, deriv)
    end

    BasisFunction1D(basis, index) = BasisFunction1D(basis, index, 0)
end

deriv(b::BasisFunction1D) = typeof(b)(b.basis, b.index, b.deriv+1)
deriv(b::BasisFunction1D, order) = typeof(b)(b.basis, b.index, b.deriv+order)

function (b::Basis1D)(pt::T) where {T<:Real}
    rng = supported(b, pt)
    (dropdims(evaluate_raw(b, [pt], b.deriv, rng), dims=2), rng)
end

function (b::BasisFunction1D)(pt::T) where {T<:Real}
    rng = supported(b.basis, pt)
    if b.index ∉ rng return 0.0 end
    evaluate_raw(b.basis, [pt], b.deriv, rng)[1 + b.index - rng.start, 1]
end

function (b::Basis1D)(pts::Vector{T}) where {T<:Real}
    tp = Tuple{Vector{T}, UnitRange{Int}}
    res = Array{tp}(undef, length(pts))

    j = 0
    for (subpts, rng) in supported(b, pts)
        out = evaluate_raw(b, subpts, b.deriv, rng)
        for i in 1:length(subpts)
            res[j+=1] = (out[:,i], rng)
        end
    end

    res
end

function (b::BasisFunction1D)(pts::Vector{T}) where {T<:Real}
    res = zeros(Float64, length(pts))

    i = 1
    for (subpts, rng) in supported(b.basis, pts)
        i += length(subpts)
        if b.index ∉ rng continue end

        out = evaluate_raw(b.basis, subpts, b.deriv, rng)
        res[i-length(subpts):i-1] = out[findall(in(b.index), rng), :]
    end

    res
end

function (b::Basis1D)(pt::S, coeffs::Vector{T}) where {S<:Real, T<:Real}
    (vals, idxs) = b(pt)
    dot(vals, coeffs[idxs])
end

function (b::Basis1D)(pt::S, coeffs::Matrix{T}) where {S<:Real, T<:Real}
    (vals, idxs) = b(pt)
    vals' * coeffs[idxs,:]
end

(b::Basis1D)(pts::Vector{S}, coeffs::Vector{T}) where {S<:Real, T<:Real} =
    Float64[dot(vals, coeffs[idxs]) for (vals, idxs) in b(pts)]

function (b::Basis1D)(pts::Vector{S}, coeffs::Matrix{T}) where {S<:Real, T<:Real}
    res = zeros(Float64, length(pts), size(coeffs, 2))

    for (i, (vals, idxs)) in enumerate(b(pts))
        res[i,:] = vals' * coeffs[idxs,:]
    end

    res
end


# Include specific bases
# ========================================================================

include("BSplineBasis.jl")
include("NURBSBasis.jl")


end  # module Bases
