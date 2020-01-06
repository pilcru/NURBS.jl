using NURBS.Bases

export Curve, BCurve, NCurve, evalpoint, evaltangent, evalderiv1

struct Curve{B<:Basis1D, T<:Real} <: Object1D
    basis::B
    d1b::B
    points::Array{T, 2}

    function Curve(basis::B, points::Array{T, 2}) where {B<:Basis1D, T<:Real}
        d1b = deriv(basis)
        new{B, T}(basis, d1b, points)
    end
end

const BCurve = Curve{BSplineBasis}
const NCurve = Curve{NURBSBasis}

function evalpoint(crv::Curve, u)
    return crv.basis(u, crv.points)
end

function evaltangent(crv::Curve, u)
    t = crv.d1b(u, crv.points)'
    t ./= norm(t)
end

function evalderiv1(crv::Curve, u)
    return crv.d1b(u, crv.points)'
end

