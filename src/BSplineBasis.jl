import IterTools: groupby, imap
using Zygote: @nograd

struct BSplineBasis <: Basis1D
    knots::Vector{Float64}
    order::Int
    deriv::Int

    function BSplineBasis(knots, order, deriv=0; extend=true)
        if deriv >= order throw(ArgumentError("Differentiation order not supported")) end

        for (kn, kp) in zip(knots[2:end], knots[1:end-1])
            if kn < kp throw(ArgumentError("Knot vector must be nondecreasing")) end
        end

        if extend
            knots = [fill(knots[1], order-1); knots; fill(knots[end], order-1)]
        else
            for d in 1:order-1
                if knots[1+d] != knots[1] || knots[end-d] != knots[end]
                    throw(ArgumentError("Expected $order repeated knots on either end"))
                end
            end
        end

        new(knots, order, deriv)
    end

    BSplineBasis(lft::Real, rgt::Real, elements::Int, order::Int) =
        BSplineBasis(range(lft, stop=rgt, length=elements+1), order)
end

const BSpline = BasisFunction1D{BSplineBasis}


Base.length(b::BSplineBasis) = length(b.knots) - b.order
Base.size(b::BSplineBasis) = (length(b),)
Base.getindex(b::BSplineBasis, i) = BSpline(b, i, b.deriv)
nderivs(b::BSplineBasis) = b.order - 1

domain(b::BSplineBasis) = Interval(b.knots[1], b.knots[end])
domain(b::BSpline) = Interval(b.basis.knots[b.index], b.basis.knots[b.index+b.basis.order])

degree(b::BSplineBasis) = b.order - b.deriv - 1
degree(b::BSpline) = b.basis.order - b.deriv - 1

deriv(b::BSplineBasis) = BSplineBasis(b.knots, b.order, b.deriv+1; extend=false)
deriv(b::BSplineBasis, order) = BSplineBasis(b.knots, b.order, b.deriv+order; extend=false)


function supported(b::BSplineBasis, pt::T) where {T<:Real}
    kidx = b.order - 1 + searchsorted(b.knots[b.order:end], pt).stop
    stop = b.knots[kidx] == b.knots[end] ? kidx - b.order : kidx
    stop - b.order + 1 : stop
end

function supported(b::BSplineBasis, pts::Vector{T}) where {T<:Real}
    (min, max) = extrema(pts)
    @assert(min in domain(b) && max in domain(b), "Domain error")

    idxs = zeros(Int, length(pts))

    if !issorted(pts)
        for (i, pt) in enumerate(pts)
            idxs[i] = supported(b, pt).stop
        end
    else
        kidx = b.order
        for (i, pt) in enumerate(pts)
            kidx += searchsorted(b.knots[kidx:end], pt).stop - 1
            idxs[i] = b.knots[kidx] == b.knots[end] ? kidx - b.order : kidx
        end
    end

    imap(groupby(i -> i[2], enumerate(idxs))) do i
        (pts[i[1][1]:i[end][1]], i[1][2] - b.order + 1 : i[1][2])
    end
end

@nograd function evaluate_raw(b::BSplineBasis, pts::Vector{T}, deriv::Int, rng::UnitRange{Int}) where {T<:Real}
    # Basis values of order 1 (piecewise constants)
    bvals = zeros(Float64, (b.order, length(pts)))
    bvals[end,:] .= 1.0
    
    p = b.order
    bi = rng.start + p

    # Order increment
    for k in 0:p-deriv-2
        bvals[p-k:end, :] ./= b.knots[bi:bi+k] - b.knots[bi-k-1:bi-1]

        for (i, kp, kn) in zip(p-k-1:p-1, b.knots[bi-k-2:bi-2], b.knots[bi:bi+k])
            bvals[i,:] .*= (pts .- kp)
            bvals[i,:] .+= bvals[i+1,:] .* (kn .- pts)
        end
        bvals[end,:] .*= (pts .- b.knots[bi-1])
    end

    # Differentiation
    for k = p-deriv-1:p-2
        bvals[p-k:end, :] ./= b.knots[bi:bi+k] - b.knots[bi-k-1:bi-1]

        bvals[1:end-1,:] = -diff(bvals, dims=1)
        bvals .*= k + 1
    end

    bvals
end
