module Surfaces

using .Bases

struct Surface{BU<:Basis1D, BV<:Basis1D, T<:Real} <: Object2D
    basisU::BU
    basisV::BV
    points::Array{T, 2}

    function Surface(basis, points)
        new(basis, points)
    end
end

const BSurface = Surface{BSplineBasis, BSplineBasis}
const NSurface = Surface{NURBSBasis, NURBSBasis}

function evalpoint(srf::Surface, u::Real, v::Real)
    return srf.basisU(u, srf.points) * srf.basisV(v, srf.points)
end

end  # module Curves
