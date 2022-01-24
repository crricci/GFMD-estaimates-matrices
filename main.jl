
using PyPlot
using Distances
using StaticArrays
using BenchmarkTools
using LinearAlgebra

include("L_GFD.jl")
include("R_DEBUG.jl")

function main()

    return 0 
end

function testLocations(N)
    """ random test locations at random """
    return locations = rand(N,2)
end


