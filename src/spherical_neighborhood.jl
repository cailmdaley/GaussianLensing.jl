"""
    SphericalNeighborSearcher(domain, locations, K)

A search method that finds `K` nearest neighbors in `domain`
`locations` according to `metric`.
"""
struct (SphericalNeighborSearcher <: AbstractNeighborSearcher)
  tree::GeoStatsBase.NearestNeighbors.BallTree
  k::Int
  locs::Vector{Int}
end

function SphericalNeighborSearcher(domain::AbstractDomain,
                                   locs::AbstractVector{Int},
                                   k::Int, radius=1.0)
  @assert 1 ≤ k ≤ length(locs) "number of neighbors must be smaller than number of data locations"
  @assert length(locs) ≤ npoints(domain) "number of data locations must be smaller than number of points"

  balltree = GeoStatsBase.NearestNeighbors.BallTree(
              coordinates(domain, locs),
              GeoStatsBase.Haversine(radius))

  SphericalNeighborSearcher(balltree, k, locs)
end

function search!(neighbors::AbstractVector{Int},
                 xₒ::AbstractVector{T},
                 searcher::SphericalNeighborSearcher) where
                 {T<:Real}

  k       = searcher.K
  inds, _ = GeoStatsBase.knn(searcher.tree, xₒ, k, true)
  locs    = view(searcher.locs, inds)

  @inbounds for i in 1:k
    neighbors[i] = locs[i]
  end
  # change to broadcast?

  k
end
