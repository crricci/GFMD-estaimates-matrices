
function f(p) 
    x,y = p
    return x^4 + y^4 + 4 * x^2 * y^2 
end


function test()
    Npt = 100

    # locations = testLocations(Npt^2) 
    locations = hcat([[x,y] for x in LinRange(-1,1,Npt) for y in LinRange(-1,1,Npt)]...)'
    locations += 0.001 * rand(Npt^2,2)

    x = locations[:,1]
    y = locations[:,2]
    Dx = x .- x'
    Dy = y .- y'
    Dist = sqrt.(Dx.^2 + Dy.^2)

    
    j = 3950   # for j in 1:size(locations,1)
    point = locations[j,:]
    
    ns = 50
    # ns = Npt^2
    idx_neighbors = closestStar(ns,j,Dist)

    idx_neighbors = fourQuadrantStar(j,Dx,Dy,Dist; nPtQuad = 3)

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

function testAll()
    Npt = 100

    # locations = testLocations(Npt^2) 
    locations = hcat([[x,y] for x in LinRange(-1,1,Npt) for y in LinRange(-1,1,Npt)]...)'
    locations += 0.001 * rand(Npt^2,2)

    x = locations[:,1]
    y = locations[:,2]
    Dx = x .- x'
    Dy = y .- y'
    Dist = sqrt.(Dx.^2 + Dy.^2)

    err∂x = zeros(Npt^2)
    err∂y = zeros(Npt^2)
    err∂x2 = zeros(Npt^2)
    err∂y2 = zeros(Npt^2)
    err∂xy = zeros(Npt^2)
    
    @showprogress for j in 1:size(locations,1)
        point = locations[j,:]
        
        ns = 50
        idx_neighbors = closestStar(ns,j,Dist)
        # idx_neighbors = fourQuadrantStar(j,Dx,Dy,Dist; nPtQuad = 3)ù

        neighbors = locations[idx_neighbors,:]
        d = Dist[idx_neighbors,j]
        h = Dx[idx_neighbors,j]
        k = Dy[idx_neighbors,j]

        u = [f(pt) for pt in eachrow(neighbors)]
        u₀ = f(point)

        Dᵤ = Dᵢ(neighbors,point,u,u₀,h,k,d)
        ∇f = ForwardDiff.gradient(f,point)
        Hf = ForwardDiff.hessian(f,point)

        err∂x[j] = abs(Dᵤ[1] - ∇f[1])/ abs(∇f[1])
        err∂y[j] = abs(Dᵤ[2] - ∇f[2])/ abs(∇f[2])
        err∂x2[j] = abs(Dᵤ[3] - Hf[1,1]) / abs(Hf[1,1])
        err∂y2[j] = abs(Dᵤ[4] - Hf[2,2]) / abs(Hf[2,2])
        err∂xy[j] = abs(Dᵤ[5] - Hf[1,2]) / abs(Hf[1,2])
    end

    nanmean(x) = mean(filter(isfinite, x))
    
    X = reshape(x,Npt,Npt)
    Y = reshape(x,Npt,Npt)'
    
    @show nanmean(err∂x)
    @show nanmean(err∂y)
    @show nanmean(err∂x2)
    @show nanmean(err∂y2)
    @show nanmean(err∂xy)

    


    return nothing 
end
