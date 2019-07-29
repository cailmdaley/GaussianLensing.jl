function covariance(Δβ::Real, Cℓs::Vector)
	cosΔβ = cos(Δβ)
	count = 0.0
	for ℓ in eachindex(Cℓs) .- 1
		count += 2ℓ / (4π) * Cℓs[ℓ+1] * legendre(cosΔβ, ℓ)
	end
	return count
end

abstract type CMBQuantity end
abstract type 𝚯 <: CMBQuantity end
abstract type E <: CMBQuantity end
abstract type B <: CMBQuantity end
abstract type ϕ <: CMBQuantity end

abstract type AbstractCMBVariogram{Q <: CMBQuantity, T <: Real,
								   D <: Metric} <: Variogram{T,D} end

struct CMBVariogram{Q, T, D} <: AbstractCMBVariogram{Q, T, D}
	Cℓs::Vector{T}
	σ₀²::T
	distance::D
end

function CMBVariogram{Q}(Cℓs::Vector{T}) where {Q, T}
	CMBVariogram{Q, T, Haversine}(Cℓs, covariance(0, Cℓs), Haversine(1.0))
end

(γ::CMBVariogram)(Δβ::Float64) = γ.σ₀² - covariance(Δβ, γ.Cℓs)
(γ::CMBVariogram)(βx, βy) = γ(evaluate(γ.distance, βx, βy))
isstationary(::CMBVariogram) = true

function cache(γ::CMBVariogram, Δβs, path)
	variogram_vals  = [γ.σ₀² - covariance(Δβ, γ.Cℓs) for Δβ in Δβs]
	writedlm(path, [Δβs variogram_vals])
end


struct CachedCMBVariogram{Q, T, D} <: AbstractCMBVariogram{Q, T ,D}
	γinterpolator::ScaledInterpolation{T}
	distance::D
end

function CachedCMBVariogram{Q}(lookup_table::Matrix{T}) where {T,Q}
	Δβs = range(0, lookup_table[end,1], length=size(lookup_table)[1])
	γs  = lookup_table[:,2]

	γinterpolator         = interpolate(γs, BSpline(Linear()))
	scaled_γinterpolator  = scale(γinterpolator, Δβs)

	return CachedCMBVariogram{Q, T, Haversine}(scaled_γinterpolator,
											   Haversine(1.0))
end

function CachedCMBVariogram{Q}(path::String) where Q <: CMBQuantity
	CachedCMBVariogram{Q}(readdlm(path))
end
(γ::CachedCMBVariogram)(Δβ::Float64) = γ.γinterpolator(Δβ)
(γ::CachedCMBVariogram)(βx, βy)      = γ(evaluate(γ.distance, βx, βy))
isstationary(::CachedCMBVariogram)   = true

@estimsolver CMBKriging begin
  @param variogram = CachedCMBVariogram{𝚯}(
    "/home/cailmdaley/spt/datafiles/CMB_variogram_lookup.txt")
  @param K = 10
end

function preprocess(problem::EstimationProblem, solver::CMBKriging)
  # retrieve problem info
  pdata = data(problem)
  pdomain = domain(problem)

  # result of preprocessing
  preproc = Dict{Symbol,NamedTuple}()

  for (var, V) in variables(problem)
    # get user parameters
    varparams = var ∈ keys(solver.params) ? solver.params[var] : CMBKrigingParam()

    # determine which Kriging variant to use
    estimator = OrdinaryKriging(varparams.variogram)

    path = SimplePath(pdomain)

    # determine maximum number of conditioning neighbors
    K = varparams.K
    if varparams.K ≠ nothing
      # locations with data for given variable
      datalocs = collect(keys(datamap(problem, var)))
	  searchdomain = domain(pdata)
      # spherical K nearest neighbor searcher
      neighsearcher = SphericalNeighborSearcher(searchdomain, datalocs, K)
    else
      # use all data points as neighbors
      neighsearcher = nothing
    end

    # save preprocessed input
    preproc[var] = (estimator=estimator, path=path,
										K=K, neighsearcher=neighsearcher)
  end

  preproc
end

function values!(buff::AbstractVector, spatialdata::AbstractSpatialData,
								 locations::AbstractVector{Int}, var::Symbol)
  for j in 1:length(locations)
		buff[j] = value(spatialdata, locations[j], var)
  end
end

function GeoStatsBase.solve(problem::EstimationProblem, solver::CMBKriging)
  # preprocess user input
  preproc = preprocess(problem, solver)

  # results for each variable
  μs = []; σs = []

  for (var, V) in variables(problem)
    if preproc[var].K ≠ nothing
      # perform Kriging with reduced number of neighbors
      varμ, varσ = solve_locally(problem, var, preproc)
    else
      # perform Kriging with all data points as neighbors
      varμ, varσ = solve_globally(problem, var, preproc)
    end

    push!(μs, var => varμ)
    push!(σs, var => varσ)
  end

  EstimationSolution(domain(problem), Dict(μs), Dict(σs))
end

function solve_locally(problem::EstimationProblem, var::Symbol, preproc)
    # retrieve problem info
    pdata = data(problem)
    pdomain = domain(problem)

    # unpack preprocessed parameters
    estimator, path, K, neighsearcher = preproc[var]

    # determine value type
    V = variables(problem)[var]

    # pre-allocate memory for result
    varμ = Vector{V}(undef, npoints(pdomain))
    varσ = Vector{V}(undef, npoints(pdomain))

    # pre-allocate memory for coordinates

    xₒ = GeoStatsBase.MVector{ndims(pdomain),coordtype(pdomain)}(undef)

    # pre-allocate memory for neighbors
    neighbors = Vector{Int}(undef, K)
    X = Matrix{coordtype(pdata)}(undef, ndims(pdata), K)
		z = Vector{V}(undef, K)

    # estimation loop
    for location in path
      coordinates!(xₒ, pdomain, location)

      # find neighbors with previously estimated values
      search!(neighbors, xₒ, neighsearcher)
      # update neighbors coordinates
      coordinates!(X, pdata, neighbors)
			values!(z, pdata, neighbors, var)

      # fit estimator to data
      krig = fit(estimator, X, z)

      # save mean and variance
      μ, σ² = predict(krig, xₒ)
      varμ[location] = μ
      varσ[location] = σ²
    end

    varμ, varσ
end

function solve_globally(problem::EstimationProblem, var::Symbol, preproc)
    # retrieve problem info
    pdata = data(problem)
    pdomain = domain(problem)

    # unpack preprocessed parameters
    estimator, path, K, neighsearcher = preproc[var]

    # determine value type
    V = variables(problem)[var]

    # pre-allocate memory for result
    varμ = Vector{V}(undef, npoints(pdomain))
    varσ = Vector{V}(undef, npoints(pdomain))

    # pre-allocate memory for coordinates
    xₒ = GeoStatsBase.MVector{ndims(pdomain),coordtype(pdomain)}(undef)

    # fit estimator to data
    X, z = valid(pdata, var)
    krig = fit(estimator, X, z)

    for location in path
      coordinates!(xₒ, pdomain, location)

      μ, σ² = predict(krig, xₒ)

      varμ[location] = μ
      varσ[location] = σ²
    end

    varμ, varσ
end

function gp_interpolate(y::Dict, X::Matrix, Xₒ::Matrix, maxneighbors=10)
	length(y) != 1 && throw("Multiple fields not implemented :(")
	var = collect(keys(y))[1]

	pdata    = PointSetData(Dict(y), X)
	pdomain  = PointSet(Xₒ)
	problem = EstimationProblem(pdata, pdomain, var,
								mapper=CopyMapper())
	solver = CMBKriging()

	return solve(problem, solver)
end
