

function allDist(locations)
    """ matrix of all pairwise distances """
    Dist = pairwise(Euclidean(),locations,dims=1)
    return Dist
end

function Dᵢ(neighbors, point, u, u₀)
    """ Compute D(point) = (∂x,∂y,∂xx,∂yy,∂xy)(point) wrt neighbors 
        u = [value at point, values of the function at the neighbors]'
    """ 

    h = - (point[1] .- neighbors[:,1])
    k = - (point[2] .- neighbors[:,2])
    d = sqrt.(h.^2 + k.^2)
    @show R = maximum(d)  # Rᵢ
    w = [weight(dᵢ,R) for dᵢ in d]

    A = compute_A(h,k,w)
    B = compute_B(h,k,w)
    b = compute_b(h,k,w,u,u₀)

    # D = inv(A) * B
    # Dᵤ = D * [u₀; u]
    # return Dᵤ
    
    return Dᵤ = A \ b
end

function compute_A(h,k,w)
    """ Create the matrix A """
    ns = length(w)  # common h,k,w
    A = zeros(5,5)

    w2 = [w[i]^2 for i in 1:ns]
    
    for i in 1:ns
        A[1,1] += h[i]^2 * w2[i]
        A[1,2] += h[i] * k[i] * w2[i]
        A[1,3] += k[i]^3 * w2[i] / 2
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
    A += triu(A,1)'     # symmetrize

    return A
end

function compute_B(h,k,w)
    """ Create the matrix B """
    ns = length(w)  # common h,k,w
    B = zeros(5,ns + 1)

    # w2 = @SVector [wᵢ^2 for wᵢ in w]
    w2 = [w[i]^2 for i in 1:ns]

    B[1,1] = - sum([h[i] * w2[i] for i in 1:ns])
    B[1,2] = - sum([k[i] * w2[i] for i in 1:ns])
    B[1,1] = - sum([h[i]^2 * w2[i] for i in 1:ns]) / 2
    B[1,1] = - sum([k[i]^2 * w2[i] for i in 1:ns]) / 2
    B[1,1] = - sum([h[i] * k[i] * w2[i] for i in 1:ns])
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

