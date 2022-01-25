
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
        idx_neighbors = closestStar(ns,j,Dist)
        # idx_neighbors = fourQuadrantStar(j,Dx,Dy,Dist; nPtQuad = 3)

        neighbors = locations[idx_neighbors,:]
        d = Dist[idx_neighbors,j]
        h = Dx[idx_neighbors,j]
        k = Dy[idx_neighbors,j]

        u = [f(pt) for pt in eachrow(neighbors)]
        u₀ = f(point)

        Du = Dᵤ(u,u₀,h,k,d)
        ∇f = ForwardDiff.gradient(f,point)
        Hf = ForwardDiff.hessian(f,point)

        err∂x[j] = abs(Du[1] - ∇f[1])/ abs(∇f[1])
        err∂y[j] = abs(Du[2] - ∇f[2])/ abs(∇f[2])
        err∂x2[j] = abs(Du[3] - Hf[1,1]) / abs(Hf[1,1])
        err∂y2[j] = abs(Du[4] - Hf[2,2]) / abs(Hf[2,2])
        err∂xy[j] = abs(Du[5] - Hf[1,2]) / abs(Hf[1,2])
    end

    nanmean(x) = mean(filter(isfinite, x))
    @show nanmean(err∂x)
    @show nanmean(err∂y)
    @show nanmean(err∂x2)
    @show nanmean(err∂y2)
    @show nanmean(err∂xy)
    
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

    Mx, My, Mxx, Myy, Mxy = compute_MDiff(locations)

    x = locations[:,1]
    y = locations[:,2]
    Dx = x .- x'
    Dy = y .- y'
    Dist = sqrt.(Dx.^2 + Dy.^2)

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

    @showprogress for i in 1:size(locations,1)
        point = locations[i,:]
        
        ns = 20
        idx_neighbors = closestStar(ns,i,Dist)

        neighbors = locations[idx_neighbors,:]
        d = Dist[idx_neighbors,i]
        h = Dx[idx_neighbors,i]
        k = Dy[idx_neighbors,i]

        u = [f(pt) for pt in eachrow(neighbors)]
        u₀ = f(point)

        Du = Dᵤ(u,u₀,h,k,d)
        # Dᵢ = compute_Dᵢ(h,k,d)
        # Du = Dᵢ * [u₀; u]

        err∂x[i] = abs(Du[1] - ∂xf[i])
        err∂y[i] = abs(Du[2] - ∂yf[i])
        err∂x2[i] = abs(Du[3] - ∂xxf[i]) 
        err∂y2[i] = abs(Du[4] - ∂yyf[i]) 
        err∂xy[i] = abs(Du[5] - ∂xyf[i]) 
    end

    nanmean(x) = mean(filter(isfinite, x))
    @show maximum(err∂x)
    @show maximum(err∂y)
    @show maximum(err∂x2)
    @show maximum(err∂y2)
    @show maximum(err∂xy)

    return nothing
end