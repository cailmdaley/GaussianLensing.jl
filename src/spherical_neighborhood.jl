"""
    SphericalNeighborSearcher(domain, locations, K)

A search method that finds `K` nearest neighbors in `domain`
`locations` according to `metric`.
"""

struct (SphericalNeighborSearcher 
        <: AbstractNeighborSearcher)
  tree::BallTree
  K::Int
  locs::Vector{Int}
end

function SphericalNeighborSearcher(domain::AbstractDomain, 
                                   locs::AbstractVector{Int}, 
                                   K::Int, radius=1.0)
  @assert 1 ≤ K ≤ length(locs) "number of neighbors must be smaller than number of data locations"
  @assert length(locs) ≤ npoints(domain) "number of data locations must be smaller than number of points"
  
  balltree = BallTree(coordinates(domain, locs), 
                      GeoStatsBase.Haversine(radius))
  SphericalNeighborSearcher(balltree, K, locs)
end

function search!(neighbors::AbstractVector{Int}, 
  
                 xₒ::AbstractVector{T},
                 searcher::SphericalNeighborSearcher) where 
                 {T<:Real}
                 
  K       = searcher.K
  inds, _ = knn(searcher.tree, xₒ, K, true)
  locs    = view(searcher.locs, inds)

  @inbounds for i in 1:K
    neighbors[i] = locs[i]
  end
  # change to broadcast?

  K
end
