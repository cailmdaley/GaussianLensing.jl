abstract type Plane end
abstract type Source <: Plane end
abstract type Image <: Plane end

abstract type State end
abstract type Unlensed  <: State end
abstract type Delensed  <: State end
abstract type   Lensed  <: State end

# struct Lens{T <: Number, P <: Plane, S <: State}
# 	m::AbstractMap{T}
# end
# Lens{P <: Plane, S<:State}(m::AbstractMap{T}) = Lens{T,P,S}(m)
# Lens{Source}(m::AbstractMap) = Lens{Source,Unlensed}(m)
# Lens{Image}(m::AbstractMap) = Lens{Image,Lensed}(m)
# 
# 
# 
# function lens(l::Lens{T,P,S}, maxneighbors=10) where {T, P<:Source, S}
# 	# lensed_data_coords = get
# 	solution = solve_kriging(:CMB => hpl.data_values, lensed_data_coords, 
# 				             hpl.domain_coords, maxneighbors) 
# 
# 	hpl.domain_mean     .= solution[:CMB].mean
# 	hpl.domain_variance .= solution[:CMB].variance
# 	return hpl.domain_mean, hpl.domain_variance
# end
