
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

const degree2km = 110   # one land degree to km 

function main()
    muniCSV = CSV.read("../coordinates_municipalities.csv",DataFrame)
    locations = hcat(muniCSV.V1,muniCSV.V2)
    (Mx2, My2, Mxx2, Myy2, Mxy2),toRemove = compute_MDiff(locations,radius = 50/degree2km,method = "dist")

    CSV.write("Mx.csv",Tables.table(vec(Mx)),header=false)
    CSV.write("My.csv",Tables.table(vec(My)),header=false)
    CSV.write("Mxx.csv",Tables.table(vec(Mxx)),header=false)
    CSV.write("Myy.csv",Tables.table(vec(Myy)),header=false)
    CSV.write("Mxy.csv",Tables.table(vec(Mxy)),header=false)
    
end

