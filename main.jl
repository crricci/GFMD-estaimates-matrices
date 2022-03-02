
using Distances
using StaticArrays
using BenchmarkTools
using LinearAlgebra
using ForwardDiff
using ProgressMeter
using CSV, DataFrames

include("L_GFD.jl")
include("L_AggregateMatrix.jl")

include("R_DEBUG.jl")

function main()
    muniCSV = CSV.read("coordinates_municipalities.csv",DataFrame)
    locations = hcat(muniCSV.V1,muniCSV.V2)
    Mx, My, Mxx, Myy, Mxy = compute_MDiff(locations)

    CSV.write("Mx.csv",Tables.table(vec(Mx)),header=false)
    CSV.write("My.csv",Tables.table(vec(My)),header=false)
    CSV.write("Mxx.csv",Tables.table(vec(Mxx)),header=false)
    CSV.write("Myy.csv",Tables.table(vec(Myy)),header=false)
    CSV.write("Mxy.csv",Tables.table(vec(Mxy)),header=false)
    
end

function testLocations(N)
    """ random test locations at random """
    return locations = rand(N,2)
end


