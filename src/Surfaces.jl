using NURBS.Bases
import LinearAlgebra: dot, cross, norm

struct Surface{BU<:Basis1D, BV<:Basis1D, T<:Real} <: Object2D
    basisU::BU
    basisV::BV
    d1bU::BU
    d2bU::BU
    d3bU::BU
    d1bV::BV
    d2bV::BV
    d3bV::BV
    points::Array{T, 3}

    function Surface(basisU::BU, basisV::BV, points::Array{T, 3}) where {BU<:Basis1D, BV<:Basis1D, T<:Real}
        d1bU = deriv(basisU)
        d2bU = deriv(d1bU)
        d3bU = deriv(d2bU)
        d1bV = deriv(basisV)
        d2bV = deriv(d1bV)
        d3bV = deriv(d2bV)
        new{BU, BV, T}(basisU, basisV, d1bU, d2bU, d3bU, d1bV, d2bV, d3bV, points)
    end
end

const BSurface = Surface{BSplineBasis, BSplineBasis}
const NSurface = Surface{NURBSBasis, NURBSBasis}

function evalpoint(srf::Surface, u, v)
    return srf.basisV(v, srf.basisU(u, srf.points))'
end

function evalnormal(srf::Surface, u, v)
    (du, dv, duu, duv, dvv) = evalderiv2(srf, u, v)
    E = dot(du, du)
    F = dot(du, dv)
    G = dot(dv, dv)
    n = cross(du, dv)
    n ./= norm(n)
    L = dot(duu, n)
    M = dot(duv, n)
    N = dot(dvv, n)
    dS2 = E*G - F^2
    K = (L*N - M^2) / dS2
    H = (E*N - 2*F*M + G*L) / (2*dS2)
    Δ = sqrt(H^2 - K)
    if H > 0  # Rhino orders principals by absolute value
        κ1 = H + Δ
        κ2 = H - Δ
    else
        κ1 = H - Δ
        κ2 = H + Δ
    end
    λ1 = (κ1*F - M) / (N - κ1*G)
    λ2 = (κ2*F - M) / (N - κ2*G)
    d1 = (du + λ1*dv)
    d1 ./= norm(d1)
    d2 = (du + λ2*dv)
    d2 ./= norm(d2)
    (n, κ1, κ2, du, dv, d1, d2, dS2)
end

function evalderiv1(srf::Surface, u, v)
    du = srf.basisV(v, srf.d1bU(u, srf.points))
    dv = srf.d1bV(v, srf.basisU(u, srf.points))
    (du', dv')
end

function evalderiv2(srf::Surface, u, v)
    duu = srf.basisV(v, srf.d2bU(u, srf.points))
    du = srf.basisV(v, srf.d1bU(u, srf.points))
    duv = srf.d1bV(v, srf.d1bU(u, srf.points))
    dv = srf.d1bV(v, srf.basisU(u, srf.points))
    dvv = srf.d2bV(v, srf.basisU(u, srf.points))
    (du', dv', duu', duv', dvv')
end

function evalderiv3(srf::Surface, u, v)
    duuu = srf.basisV(v, srf.d3bU(u, srf.points))
    duu = srf.basisV(v, srf.d2bU(u, srf.points))
    duuv = srf.d1bV(v, srf.d2bU(u, srf.points))
    duv = srf.d1bV(v, srf.d1bU(u, srf.points))
    du = srf.basisV(v, srf.d1bU(u, srf.points))
    duvv = srf.d2bV(v, srf.d1bU(u, srf.points))
    dv = srf.d1bV(v, srf.basisU(u, srf.points))
    dvv = srf.d2bV(v, srf.basisU(u, srf.points))
    dvvv = srf.d3bV(v, srf.basisU(u, srf.points))
    (du', dv', duu', duv', dvv', duuu', duuv', duvv', dvvv')
end

function BSurface(orderU::Int, orderV::Int, points::Array{T, 3}) where {T<:Real}
    bu = BSplineBasis(range(-1, stop=1, length=size(points, 1)-orderU+2), orderU)
    bv = BSplineBasis(range(-1, stop=1, length=size(points, 2)-orderV+2), orderV)
    Surface(bu, bv, points)
end

function BSurface(KnotsU::Vector{T}, KnotsV::Vector{T}, orderU::Int, orderV::Int, points) where {T<:Real}
    bu = BSplineBasis(2 .* (KnotsU .- minimum(KnotsU)) ./ (maximum(KnotsU) - minimum(KnotsU)) .- 1, orderU)
    bv = BSplineBasis(2 .* (KnotsV .- minimum(KnotsV)) ./ (maximum(KnotsV) - minimum(KnotsV)) .- 1, orderV)
    Surface(bu, bv, points)
end

function isocurve(srf::Surface, dir::Bool, u)
    b = dir ? srf.basisV : srf.basisU
    p = dir ? srf.basisU(u, srf.points) : srf.basisV(u, permutedims(srf.points, [2, 1, 3]))
    return Curve(b, p)
end

NSurface(srf::BSurface) = Surface(NURBSBasis(srf.basisU), NURBSBasis(srf.basisV))
NSurface(points::Array{T, 3}, degU::Int, degV::Int) where {T<:Real} = NSurface(BSurface(points, degU, degV))
