module GaussianLensing

using Reexport 
@reexport using Unitful, UnitfulAstro
@reexport using Plots, LaTeXStrings

############################################################################
# Utilties & Pycall

using DelimitedFiles
include("utils.jl")
export read_Cℓs

@reexport using Healpix
# include("maps.jl")
# export 
# 	MaskedMap,
# 	makeCMB, makeCMBmasked, 
# 	make_αs, make_αs_masked,
# 	reorder, 
# 	getcoords,
# 	plotCMB
	
# using PyCall; const healpy = PyNULL()
# function __init__()
#     copy!(healpy, pyimport("healpy"))
# end
# export healpy

############################################################################
# Kriging & GeoStats 
@reexport using GeoStats
using GeoStatsBase 

using NearestNeighbors
include("spherical_neighborhood.jl")

using Jacobi
using Interpolations
include("kriging.jl")
export 
	# CMBVariogram,
	solve_kriging

	
############################################################################
# Lensing

# include("lensing.jl")
# export
# 	HealpixLens,
# 	lens

	
end # module GaussianLensing
