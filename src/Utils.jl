module Utils

export Interval

mutable struct Interval{T<:Real}
    lower::T
    upper::T
end

Base.in(x::T, interval::Interval{T}) where {T<:Real} = interval.lower <= x <= interval.upper
Base.:(==)(a::Interval, b::Interval) = a.lower == b.lower && a.upper == b.upper
⊆(a::Interval, b::Interval) = b.lower <= a.lower <= a.upper <= b.upper
⊈(a::Interval, b::Interval) = !(a ⊆ b)
⊊(a::Interval, b::Interval) = a ⊆ b && (b.lower < a.lower || a.upper < b.upper)

end  # module Utils
