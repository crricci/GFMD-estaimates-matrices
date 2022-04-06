
function f(p) 
    x,y = p
    return x^4 + y^4 + 4 * x^2 * y^2 
end

function testAllErrors()
    Npt = 100

    locations = testLocations(Npt^2) 
    # locations = hcat([[x,y] for x in LinRange(-1,1,Npt) for y in LinRange(-1,1,Npt)]...)'
    # locations += 0.001 * rand(Npt^2,2)

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
        
        ns = 20
        radius = 0.25
        idx_neighbors = closestStar(ns,j,Dist)
        # idx_neighbors = fourQuadrantStar(j,Dx,Dy,Dist; nPtQuad = 3)
        # idx_neighbors = radiusStar(radius,j,Dist)   # locations within radius

        neighbors = locations[idx_neighbors,:]
        d = Dist[idx_neighbors,j]
        h = Dx[idx_neighbors,j]
        k = Dy[idx_neighbors,j]

        u = [f(pt) for pt in eachrow(neighbors)]
        u₀ = f(point)

        Du = Dᵤ(u,u₀,h,k,d)
        ∇f = ForwardDiff.gradient(f,point)
        Hf = ForwardDiff.hessian(f,point)

        err∂x[j] = abs(Du[1] - ∇f[1])
        err∂y[j] = abs(Du[2] - ∇f[2])
        err∂x2[j] = abs(Du[3] - Hf[1,1]) 
        err∂y2[j] = abs(Du[4] - Hf[2,2]) 
        err∂xy[j] = abs(Du[5] - Hf[1,2]) 
    end

    nanmean(x) = mean(filter(isfinite, x))
    @show nanmean(err∂x)
    @show nanmean(err∂y)
    @show nanmean(err∂x2)
    @show nanmean(err∂y2)
    @show nanmean(err∂xy)
    

    # x = locations[:,1]
    # y = locations[:,2]
    # Dx = x .- x'
    # Dy = y .- y'
    # Dist = sqrt.(Dx.^2 + Dy.^2)
    # X = reshape(x,Npt,Npt)
    # Y = reshape(x,Npt,Npt)'
    # surf(X,Y,reshape(err∂y,Npt,Npt),cmap=ColorMap("summer"))

    return nothing 
    
end

function testAllAggregateDiff()
    Npt = 100
    
    locations = testLocations(Npt^2) 
    # locations = hcat([[x,y] for x in LinRange(-1,1,Npt) for y in LinRange(-1,1,Npt)]...)'
    # locations += 0.001 * rand(Npt^2,2)

    (Mx, My, Mxx, Myy, Mxy),toRemove = compute_MDiff(locations; method = "dist", radius = 0.25)

    u = [f(p) for p in eachrow(locations)]
    ∂xf = Mx * u
    ∂yf = My * u
    ∂xxf = Mxx * u
    ∂yyf = Myy * u
    ∂xyf = Mxy * u
    
    err∂x = zeros(Npt^2)
    err∂y = zeros(Npt^2)
    err∂x2 = zeros(Npt^2)
    err∂y2 = zeros(Npt^2)
    err∂xy = zeros(Npt^2)

    for i in 1:size(locations,1)
        point = locations[i,:]
        
        ∇f = ForwardDiff.gradient(f,point)
        Hf = ForwardDiff.hessian(f,point)

        err∂x[i] = abs(∂xf[i] - ∇f[1])
        err∂y[i] = abs(∂yf[i] - ∇f[2])
        err∂x2[i] = abs(∂xxf[i] - Hf[1,1]) 
        err∂y2[i] = abs(∂yyf[i] - Hf[2,2]) 
        err∂xy[i] = abs(∂xyf[i] - Hf[1,2]) 
    end

    nanmean(x) = mean(filter(isfinite, x))
    @show nanmean(err∂x)
    @show nanmean(err∂y)
    @show nanmean(err∂x2)
    @show nanmean(err∂y2)
    @show nanmean(err∂xy)

    # x = locations[:,1]
    # y = locations[:,2]
    # Dx = x .- x'
    # Dy = y .- y'
    # Dist = sqrt.(Dx.^2 + Dy.^2)
    # X = reshape(x,Npt,Npt);
    # Y = reshape(x,Npt,Npt)';
    # surf(X,Y,reshape(err∂y,Npt,Npt),cmap=ColorMap("summer"))

    return nothing
end


function testLocations(N)
    """ random test locations at random """
    return locations = rand(N,2)
end

