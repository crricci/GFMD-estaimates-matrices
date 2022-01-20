
function testLocations(N)
    """random test locations at random
    """
    return locations = rand(N,2)
end

function allDist(locations)
    x = @view locations[:,1]
    y = @view locations[:,2]
    D = pairwise(Euclidean(),x,y)
    return D
end

function Dᵢ(neighbors, point)
    """ Compute D(point) = (∂x,∂y,∂xx,∂yy,∂xy)(point) wrt neighbors
    """
    ns = size(neighbors,2)
    h = point[1] .- neighbors[:,1]
    k = point[2] .- neighbors[:,2]
    d = sqrt(h.^2 + y.^2)
    w 
end

function w(d,dₘ)
    x = d/dₘ
    if x ≤ 1
        return 1 - 6 * x^2 + 8 * x^3 - 3 * x^4
    else 
        return 0
    end
end

