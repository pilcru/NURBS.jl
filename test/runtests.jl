using Test
using NURBS

@testset "NURBS.jl" begin
    @testset "BSplineBasis.jl" begin
        include("BSplineBasis.jl")
    end
    @testset "NURBSBasis.jl" begin
        include("NURBSBasis.jl")
    end
end
