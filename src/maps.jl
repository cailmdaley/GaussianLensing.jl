###########################################################################
# Map & MaskedMap initialization functions

function get_latlong_deg(m::OrderedMap, ind) 
	colat, lon = pix2ang(m, ind)
	rad2deg.((colat2lat(colat), lon))
end

struct MaskedMap{T, O} <: OrderedMap{T, O}
	pixels::Vector{T}
	inds::Vector{Int64}
	resolution::Resolution
end

# Map to MaskedMap
function MaskedMap(m::Map{T, O}, θlims, ϕlims) where {T,O}
	inds = Vector{Int64}()
	for i in eachindex(m.pixels)
		θ, ϕ  = get_latlong_deg(m, i)
		if (θlims[1] <= θ <= θlims[2]) & (ϕlims[1] <= ϕ <= ϕlims[2])
			push!(inds, i)
		end
	end
	MaskedMap{T,O}(m.pixels[inds], inds, m.resolution)
end

# MaskedMap to Map
function Healpix.Map(mm::MaskedMap{T, O}) where {T, O}
	healpixels = fill(NaN, mm.resolution.numOfPixels)
	healpixels[mm.inds] = mm.pixels
	Map{O}(healpixels)
end

# Pixel vector to MaskedMap
function MaskedMap{O}(healpixels::Vector{T}, θlims, ϕlims) where 
					  {T <: Number, O <: Order}
	MaskedMap(Map{O}(healpixels), θlims, ϕlims)
end

# Pixel vector to defaults to ringOrdering
function MaskedMap(healpixels::Vector{T}, θlims, ϕlims) where T <: Number
	MaskedMap{RingOrder}(healpixels, θlims, ϕlims)
end
###########################################################################
# Methods for CMB Map Making

makeCMB(nside) = Map(healpy.synfast(read_Cℓs(), nside))
makeCMB(nside, θlims, ϕlims) = MaskedMap(makeCMB(nside), θlims, ϕlims)

function make_αs(nside) 
	healpy.alm2map_der1(healpy.synalm(read_Cℓs()), nside)[2:end, :]
end
function make_αs(m::Map{T,O}) where {T, O <: RingOrder}
	make_αs(m.resolution.nside)
end
function make_αs(mm::MaskedMap{T,O}) where {T, O <: RingOrder}
	make_αs(mm.resolution.nside)[:,mm.inds]
end

# function make_αs(mm:MaskedMap{T,O}, coeff) where {T, O <: RingOrder}
# 	coeff .* healpy.alm2map_der1(healpy.synalm(
# 		read_Cℓs()), m.resolution.nside)[2:end, mm.inds]
	# end

###########################################################################
# Methods for MaskedMap reordering and coordinate querying

get_reorder_func(mm::MaskedMap{T, RingOrder}) where T = ring2nest
get_reorder_func(mm::MaskedMap{T, NestedOrder}) where T = nest2ring
function reorder(mm::MaskedMap{T,O}) where {T, O}
	new_inds = [get_reorder_func(mm)(mm.resolution, ind) for ind in mm.inds]
	sortperm!(new_inds, new_inds)
	MaskedMap{T,NestedOrder}(mm.resolution, new_inds, mm.pixels[new_inds])
end

function resize(m::Map{T,O}, nside) where {T, O <: RingOrder}
	Map(healpy.ud_grade(m.pixels, nside))
end
function resize(mm::MaskedMap{T,O}, nside) where {T, O <: RingOrder}
	θlims, ϕlims = extrema(getcoords(mm), dims=2)
	MaskedMap(resize(Map(mm), nside), θlims, ϕlims)
end

get_healpixel_inds(m::Map) = eachindex(m.pixels)
get_healpixel_inds(mm::MaskedMap) = [mm.inds[i] for i in eachindex(mm.inds)]
function getcoords(m::OrderedMap)
	coords = zeros((2, length(m)))
	healpixel_inds = get_healpixel_inds(m)
	for i in eachindex(healpixel_inds)
		coords[:,i] .= get_latlong_deg(m, healpixel_inds[i])
	end
	coords
end
###########################################################################
# Methods for plotting Masked Maps

plotCMB(CMB::Map, c=:balance) = plot(CMB, c=c)
function GaussianLensing.plotCMB(mm::MaskedMap; rot=0, center=nothing, 
								 fov_rad=nothing, c=:balance)
	coords = getcoords(mm) 
	(θmin, θmax), (ϕmin, ϕmax) = [extrema(coords[i,:]) for i in 1:2]
			   
	if center === nothing
		center = deg2rad.([-(θmin + θmax) / 2, (ϕmin + ϕmax) / 2, rot]) 
		if isapprox(θmin, -90, atol=1)
			center[1] = deg2rad(θmin); fov_rad *= 2
		elseif isapprox(θmax, 90, atol=1)
			center[1] = deg2rad(θmax); fov_rad *= 2
		end
	end
	if fov_rad === nothing
		distance = GeoStatsBase.Haversine(1.)
		fov_rad = maximum([
			evaluate(distance, (ϕmin, θmin), (ϕmin, θmax)),
			evaluate(distance, (ϕmin, θmin), (ϕmax, θmin))]) / 2
	end
		
	plot(Map(mm), gnomonic, Dict(:center => center, :fov_rad => fov_rad),
		 c=c)
end
