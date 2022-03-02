
function compute_MDiff(x,y; ns = 20)
    Mx, My, Mxx, Myy, Mxy = compute_MDiff([x y]; ns = ns)
    return [Mx,My,Mxx,Myy,Mxy]
end

function compute_MDiff(locations; ns = 20)
    """ Mx⋅f ≈ ∂x f
        My⋅f ≈ ∂y f
        Mxx⋅f ≈ ∂xx f
        Myy⋅f ≈ ∂yy f
        Mxy⋅f ≈ ∂xy f
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

    @showprogress for i in 1:NPt

        idx_neighbors = closestStar(ns,i,Dist)
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

    return  Mx, My, Mxx, Myy, Mxy
end

