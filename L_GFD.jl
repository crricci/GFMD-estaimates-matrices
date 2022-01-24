
function testLocations(N)
    """ random test locations at random """
    return locations = rand(N,2)
end

function allDist(locations)
    """ matrix of all pairwise distances """
    Dist = pairwise(Euclidean(),locations,dims=1)
    return Dist
end

function Dᵢ(neighbors, point)
    """ Compute D(point) = (∂x,∂y,∂xx,∂yy,∂xy)(point) wrt neighbors """
    ns = size(neighbors,1)
    h = point[1] .- neighbors[:,1]
    k = point[2] .- neighbors[:,2]
    d = sqrt.(h.^2 + k.^2)
    R = maximum(d)  # Rᵢ
    w = [weight(dᵢ,R) for dᵢ in d]

end

function computeA(h,k,w)
    """ Create the matrix A """
    ns = length(w)  # common h,k,w
    A = zeros(5,5)

    w2 = @SVector [wᵢ^2 for wᵢ in w]
    
    for i in 1:ns
        A[1,1] += h[i]^2 * w2[i]
        A[1,2] += h[i] * k[i] * w2[i]
        A[1,3] += k[i]^3 * w2[i] / 2
        A[1,4] += h[i] * k[i]^2 * w2[i] / 2
        A[1,5] += h[i]^2 * k[i] * w2[i] / 2
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

function computeB(h,k,w)
    """ Create the matrix B """
    ns = length(w)  # common h,k,w
    B = zeros(5,ns + 1)

    w2 = @SVector [wᵢ^2 for wᵢ in w]

    B[1,1] = - sum(@SVector [h[i] * w2[i] for i in 1:ns])
    B[1,2] = - sum(@SVector [k[i] * w2[i] for i in 1:ns])
    B[1,1] = - sum(@SVector [h[i]^2 * w2[i] for i in 1:ns]) / 2
    B[1,1] = - sum(@SVector [k[i]^2 * w2[i] for i in 1:ns]) / 2
    B[1,1] = - sum(@SVector [h[i] * k[i] * w2[i] for i in 1:ns])
    for i in 1:ns
        B[1,i+1] = h[i] * w2[i]
        B[2,i+1] = k[i] * w2[i]
        B[3,i+1] = h[i]^2 * w2[i] / 2
        B[4,i+1] = k[i]^2 * w2[i] / 2
        B[5,i+1] = h[i] * k[i] * w2[i] 
    end

    return B
end


function weight(d,dₘ)
    x = d/dₘ
    if x ≤ 1
        return 1 - 6 * x^2 + 8 * x^3 - 3 * x^4
    else 
        return 0
    end
end

