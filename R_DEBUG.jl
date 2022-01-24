
function f(p) 
    x,y = p
    return sin(x)*cos(y)
end


function test()
    locations = testLocations(10_000) 
    Dist = allDist(locations)   # suboptimal since uses h and k but ok
    
    ns = 12 
    
    j = 40   # for j in 1:size(locations,1)
    point = locations[j,:]
    d = Dist[j,:]
    idx_neighbors = partialsortperm(d,1:ns + 1) # +1 counting itself
    
    neighbors = locations[idx_neighbors[2:end],:]
    u = [f(pt) for pt in eachrow(neighbors)]
    u = [f(point); u]


    @show Dᵤ = Dᵢ(neighbors,point,u)
    @show ForwardDiff.gradient(f,point)
    Dᵤ[1:2] ≈ ForwardDiff.gradient(f,point)
end
