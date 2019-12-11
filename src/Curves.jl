module Curves

using .Bases

struct Curve{B<:Basis1D, T<:Real} <: Object1D
    basis::B
    points::Array{T, 2}

    function Curve(basis, points)
        new(basis, points)
    end

    function Curve(points, degU, degV)
        b = BSplineBasis()
end

const BCurve = Curve{BSplineBasis}
const NCurve = Curve{NURBSBasis}

function evalpoint(crv::Curve, u::Real)
    return crv.basis(u, crv.points)
end

end  # module Curves
