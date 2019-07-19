struct (CMBVariogram{T <: Real, D <: Metric}
		<: Variogram{T,D})
	Cℓs::Vector{T}
	distance::D
	γinterpolator::Interpolations.ScaledInterpolation
end
function CMBVariogram()
	Cℓs = read_Cℓs()
	
	lookup_table = readdlm("datafiles/CMB_variogram_lookup.txt")
	γinterpolator = interpolate(lookup_table[:,2], BSpline(Linear())) 
	Δβs = range(0, lookup_table[end,1], length=size(lookup_table)[1])
	scaled_γinterpolator  = Interpolations.scale(γinterpolator, Δβs)
	
	return CMBVariogram(Cℓs, GeoStatsBase.Haversine(1.0), scaled_γinterpolator)
end
(γ::CMBVariogram)(Δβ::Float64) = γ.γinterpolator(Δβ) # γ.σ₀² - CMBcov(Δβ, γ.Cℓs)
(γ::CMBVariogram)(βx, βy) = γ(evaluate(γ.distance, βx, βy))
isstationary(::CMBVariogram)   = true

function CMBcov(Δβ::Real, Cℓs::Vector) 
	cosΔβ = cos(Δβ); count = 0.0
	for ℓ in eachindex(Cℓs)
		count += 2(ℓ-1) / (4π) * Cℓs[ℓ] * legendre(cosΔβ, ℓ-1) 
	end
	return count
end

@estimsolver CMBKriging begin
  @param variogram = CMBVariogram()
  @param K = 10
end

using GaussianLensing: SphericalNeighborSearcher

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
    
    xₒ = GeoStats.MVector{ndims(pdomain),coordtype(pdomain)}(undef)

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
      krig = KrigingEstimators.fit(estimator, X, z)

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
    xₒ = MVector{ndims(pdomain),coordtype(pdomain)}(undef)

    # fit estimator to data
    X, z = valid(pdata, var)
    krig = KrigingEstimators.fit(estimator, X, z)

    for location in path
      coordinates!(xₒ, pdomain, location)

      μ, σ² = predict(krig, xₒ)

      varμ[location] = μ
      varσ[location] = σ²
    end

    varμ, varσ
end

function solve_kriging(y, X, Xₒ, maxneighbors=10)
	pdata    = PointSetData(Dict(y), X)
	pdomain  = PointSet(Xₒ)
	var= y.first
	
	problem = EstimationProblem(pdata, pdomain, var, mapper=CoppyMapper())
	solver = CMBKriging()
	
	return solve(problem, solver)
end
