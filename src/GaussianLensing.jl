module GaussianLensing
using Reexport

############################################################################
# Utilties & Pycall

@reexport using DelimitedFiles
include("utils.jl")
export read_Cℓs

@reexport using Healpix
@reexport using Plots
include("maps.jl")
export
	MaskedMap,
	makeCMB, makeCMBmasked,
	make_αs, make_αs_masked,
	plotmap

############################################################################
# Kriging & GeoStats

@reexport using Unitful, UnitfulAstro
using GeoStatsBase
using KrigingEstimators: OrdinaryKriging, fit, predict
using Variography: Variogram
using Distances: Metric, Haversine
using Interpolations: interpolate, Linear, BSpline, ScaledInterpolation, scale
using Jacobi: legendre

include("spherical_neighborhood.jl")
include("interpolation.jl")
export
	𝚯, E, B, ϕ,
	CMBVariogram, CachedCMBVariogram, cache,
	CMBKriging, gp_interpolate



############################################################################
# Lensing

# include("lensing.jl")
# export
# 	HealpixLens,
# 	lens


end # module GaussianLensing
