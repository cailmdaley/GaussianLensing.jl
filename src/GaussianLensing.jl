module GaussianLensing
using Reexport

############################################################################
# Utilties & Pycall

using DelimitedFiles
include("utils.jl")
export read_Cℓs

@reexport using Healpix
using Plots
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
using Variography: Variogram
using Distances: Metric, Haversine
using Interpolations: interpolate, Linear, BSpline, ScaledInterpolation, scale
using Jacobi: legendre

include("spherical_neighborhood.jl")
include("kriging.jl")
export
	CMBVariogram, CachedCMBVariogram, cache,
	𝚯, E, B, ϕ,
	solve_kriging


############################################################################
# Lensing

# include("lensing.jl")
# export
# 	HealpixLens,
# 	lens


end # module GaussianLensing
