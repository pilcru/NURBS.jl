module Objects

export Object, Object1D, Object2D,
       Curve, BCurve, NCurve,
       Surface, BSurface, NSurface


# Abstract Basis types
# ========================================================================

abstract type Object end
abstract type Object1D <: Object end
abstract type Object2D <: Object end

mutable struct ObjectFunction1D{O<:Object1D}
    object::O
    index::Int
    deriv::Int

    function ObjectFunction1D{O}(object, index, deriv) where {O<:Object1D}
        if !(1 <= index <= length(object.basis)) throw(BoundsError()) end
        if !(0 <= deriv <= nderivs(object.basis)) throw(ArgumentError("Differentiation order not supported")) end
        new(object, index, deriv)
    end

    ObjectFunction1D(object, index) = ObjectFunction1D(object, index, 0)
end

mutable struct ObjectFunction2D{O<:Object2D}
    object::O
    indexU::Int
    indexV::Int
    deriv::Int

    function ObjectFunction2D{O}(object, indexU, indexV, deriv) where {O<:Object1D}
        if !(1 <= indexU <= size(object.basis, 1)) throw(BoundsError()) end
        if !(1 <= indexU <= size(object.basis, 2)) throw(BoundsError()) end
        if !(0 <= deriv <= nderivs(object.basis)) throw(ArgumentError("Differentiation order not supported")) end
        new(object, indexU, indexV, deriv)
    end

    ObjectFunction2D(object, indexU, indexV) = ObjectFunction2D(object, indexU, indexV, 0)
end

include("Curves.jl")
include("Surfaces.jl")

end  # module Objects