

function allDist(locations)
    """ matrix of all pairwise distances """
    Dist = pairwise(Euclidean(),locations,dims=1)
    return Dist
end

function Dᵤ(u, u₀, h, k, d)
    """ Compute Dᵤ(point) = (∂x,∂y,∂xx,∂yy,∂xy)(point) wrt neighbors 
        u = [value at point, values of the function at the neighbors]'
    """ 

    R = maximum(d)  # Rᵢ
    w = [weight(dᵢ,R) for dᵢ in d]

    A = compute_A(h,k,w)
    b = compute_b(h,k,w,u,u₀)

    return Dᵤ = A \ b
end

function compute_Dᵢ(h, k, d)
    """ Compute Dᵢ = matrix such that, if 
        u = [value at point, values of the function at the neighbors]'
        then Dᵤ(point) = (∂x,∂y,∂xx,∂yy,∂xy)(point) = Dᵢ * u
    """ 

    R = maximum(d)  # Rᵢ
    w = [weight(dᵢ,R) for dᵢ in d]

    A = compute_A(h,k,w)
    B = compute_B(h,k,w)

    if det(A) == 0
       # @infiltrate
    end
    
    Dᵢ = inv(A) * B
    return Dᵢ
end

function compute_A(h,k,w)
    """ Create the matrix A """
    ns = length(w)  # common h,k,w
    A = zeros(5,5)

    w2 = [w[i]^2 for i in 1:ns]
    
    for i in 1:ns
        A[1,1] += h[i]^2 * w2[i]
        A[1,2] += h[i] * k[i] * w2[i]
        A[1,3] += h[i]^3 * w2[i] / 2
        A[1,4] += h[i] * k[i]^2 * w2[i] / 2
        A[1,5] += h[i]^2 * k[i] * w2[i]
        A[2,2] += k[i]^2 * w2[i]
        A[2,3] += h[i]^2 * k[i] * w2[i] / 2
        A[2,4] += k[i]^3 * w2[i] / 2
        A[2,5] += h[i] * k[i]^2 * w2[i] 
        A[3,3] += h[i]^4 * w2[i] / 4
        A[3,4] += h[i]^2 * k[i]^2 * w2[i] / 4
        A[3,5] += h[i]^3 * k[i] * w2[i] / 2
        A[4,4] += k[i]^4 * w2[i] / 4
        A[4,5] += h[i] * k[i]^3  * w2[i] / 2
        A[5,5] += h[i]^2 * k[i]^2 * w2[i] 
    end
    A += triu(A,1)'     # symmetrize off diagonal

    return A
end

function compute_B(h,k,w)
    """ Create the matrix B """
    ns = length(w)  # common h,k,w
    B = zeros(5,ns + 1)

    w2 = [w[i]^2 for i in 1:ns]

    B[1,1] = - sum([h[i] * w2[i] for i in 1:ns])
    B[2,1] = - sum([k[i] * w2[i] for i in 1:ns])
    B[3,1] = - sum([h[i]^2 * w2[i] for i in 1:ns]) / 2
    B[4,1] = - sum([k[i]^2 * w2[i] for i in 1:ns]) / 2
    B[5,1] = - sum([h[i] * k[i] * w2[i] for i in 1:ns])
    for i in 1:ns
        B[1,i+1] = h[i] * w2[i]
        B[2,i+1] = k[i] * w2[i]
        B[3,i+1] = h[i]^2 * w2[i] / 2
        B[4,i+1] = k[i]^2 * w2[i] / 2
        B[5,i+1] = h[i] * k[i] * w2[i] 
    end

    return B
end

function compute_b(h,k,w,u,u₀)
    """ compute vector b (to remove) """ 
    ns = length(w)  # common h,k,w
    w2 = [w[i]^2 for i in 1:ns]
    
    b = zeros(5)
    for i in 1:ns
        b[1] += (- u₀ + u[i]) * h[i] * w2[i]
        b[2] += (- u₀ + u[i]) * k[i] * w2[i]
        b[3] += (- u₀ + u[i]) * h[i]^2 * w2[i] / 2
        b[4] += (- u₀ + u[i]) * k[i]^2 * w2[i] / 2
        b[5] += (- u₀ + u[i]) * h[i] * k[i] * w2[i]
    end

    return b
end


function weight(d,dₘ)
    x = d/dₘ
    if x ≤ 1
        return 1 - 6 * x^2 + 8 * x^3 - 3 * x^4
    else 
        return 0
    end
end

function closestStar(ns, j, Dist)
    """ returns the index inside locations of the ns neighbors 
        chosen by the closest point criterion, that is returns the indices
        of the ns closest locations to point
        j itself is not in idx_neighbors
    """
    d = Dist[:,j]
    idx_neighbors = partialsortperm(d,1:ns+1)[2:end] # +1 counting itself
    return idx_neighbors
end

function radiusStar(radius, j, Dist)
    """ returns the index inside locations of all the locations with distance
        smaller than "radius" from location j 
        j itself is not in idx_neighbors
    """
    d = Dist[:,j]
    idx_neighbors = findall(x -> x <= radius && x != 0,d)
    return idx_neighbors
end

function fourQuadrantStar(j,Dx,Dy,Dist; nPtQuad = 2)
    """ returns the index inside locations of the start 
        selected by the four quadrant criterion, that is returns the indices
        of the nPtQuad closes point in each of the four coordinate quadrants
        j itself is not in idx_neighbors
    """
    dx = @view Dx[:,j]
    dy = @view Dy[:,j]

    points = [(dx[i],dy[i]) for i in eachindex(dx)]
    
    I = findall(el ->  el[1] ≥ 0 && el[2] > 0, points)      # indices in I quadrant
    II = findall(el ->  el[1] < 0 && el[2] ≥ 0, points)     # indices in II quadrant
    III = findall(el ->  el[1] ≤ 0 && el[2] < 0, points)    # indices in III quadrant
    IV = findall(el ->  el[1] > 0 && el[2] ≤ 0, points)     # indices in IV quadrant
    
    idx_neighbors = Int[]
    push!(idx_neighbors, I[partialsortperm(Dist[I,j],1:nPtQuad)]...)
    push!(idx_neighbors, II[partialsortperm(Dist[II,j],1:nPtQuad)]...)
    push!(idx_neighbors, III[partialsortperm(Dist[III,j],1:nPtQuad)]...)
    push!(idx_neighbors, IV[partialsortperm(Dist[IV,j],1:nPtQuad)]...)

    return idx_neighbors
end
