
function f(p) 
    x,y = p
    return sin(x)*cos(y)
end


function test()
    Npt = 100

    locations = testLocations(Npt^2) 
    locations = hcat([[x,y] for x in LinRange(0,1,Npt) for y in LinRange(0,1,Npt)]...)'

    x = locations[:,1]
    y = locations[:,2]
    Dx = x .- x'
    Dy = y .- y'
    Dist = sqrt.(Dx.^2 + Dy.^2)

    # ns = Npt^2
    ns = 50
    
    j = 3950   # for j in 1:size(locations,1)
    point = locations[j,:]

    idx_neighbors = closestStar(ns,j,Dist)
    
    neighbors = locations[idx_neighbors,:]
    d = Dist[idx_neighbors,j]
    h = Dx[idx_neighbors,j]
    k = Dy[idx_neighbors,j]

    u = [f(pt) for pt in eachrow(neighbors)]
    u₀ = f(point)

    @show Dᵤ = Dᵢ(neighbors,point,u,u₀,h,k,d)
    @show ForwardDiff.gradient(f,point)
    @show ForwardDiff.hessian(f,point)
    
    scatter(neighbors[:,1],neighbors[:,2])
    scatter(point...)
end
