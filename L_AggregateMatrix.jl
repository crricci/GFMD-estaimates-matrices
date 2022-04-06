
function compute_MDiff(x,y;  method = "kneigh", ns = 20, radius = 50)
    (Mx, My, Mxx, Myy, Mxy), toRemove = compute_MDiff([x y]; method = method, ns =  ns, radius = radius)
    return (Mx,My,Mxx,Myy,Mxy), toRemove
end

function compute_MDiff(locations; method = "kneigh", ns = 20, radius = 50)
    """ Mx⋅f ≈ ∂x f
        My⋅f ≈ ∂y f
        Mxx⋅f ≈ ∂xx f
        Myy⋅f ≈ ∂yy f
        Mxy⋅f ≈ ∂xy f
        methods are:    "kneigh" for ns nearest neighbors 
                        "dist" for neighbors within distance radius
    """ 
    NPt = size(locations,1)

    x = locations[:,1]
    y = locations[:,2]
    Dx = x .- x'
    Dy = y .- y'
    Dist = sqrt.(Dx.^2 + Dy.^2)

    Mx = zeros(NPt,NPt) 
    My = zeros(NPt,NPt) 
    Mxx = zeros(NPt,NPt) 
    Myy = zeros(NPt,NPt) 
    Mxy = zeros(NPt,NPt) 
    
    toRemove = []
    @showprogress for i in 1:NPt
        
        if method == "kneigh" 
            idx_neighbors = closestStar(ns,i,Dist)    # ns closest locations 
        elseif method == "dist"
            idx_neighbors = radiusStar(radius,i,Dist)   # locations within radius
        else
            error("method can only be kneigh or dist")
        end 

        if length(idx_neighbors) < 6
            push!(toRemove,i)
        else 
            d = Dist[idx_neighbors,i]
            h = Dx[idx_neighbors,i]
            k = Dy[idx_neighbors,i]

            Dᵢ = compute_Dᵢ(h,k,d)
            
            Mx[i,[i;idx_neighbors]]  = Dᵢ[1,:]
            My[i,[i;idx_neighbors]]  = Dᵢ[2,:]
            Mxx[i,[i;idx_neighbors]] = Dᵢ[3,:]
            Myy[i,[i;idx_neighbors]] = Dᵢ[4,:]
            Mxy[i,[i;idx_neighbors]] = Dᵢ[5,:]
        end
    end

    Mx = Mx[setdiff(1:end, toRemove), setdiff(1:end, toRemove)]
    My = My[setdiff(1:end, toRemove), setdiff(1:end, toRemove)]
    Mxx = Mxx[setdiff(1:end, toRemove), setdiff(1:end, toRemove)]
    Myy = Myy[setdiff(1:end, toRemove), setdiff(1:end, toRemove)]
    Mxy = Mxy[setdiff(1:end, toRemove), setdiff(1:end, toRemove)]

    return  (Mx, My, Mxx, Myy, Mxy), toRemove
end

