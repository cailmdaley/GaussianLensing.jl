###########################################################################
# Methods for CMB Map Making

makeCMB(nside, Cℓs) = reorder(Map{Ring}(healpy.synfast(Cℓs, nside)))
makeCMB(nside, Cℓs, θlims, ϕlims) = MaskedMap(makeCMB(nside, Cℓs),
											  θlims, ϕlims)

function make_αs(nside)
	healpy.alm2map_der1(healpy.synalm(read_Cℓs()), nside)[2:end, :]
end
function make_αs(m::Map{T,O}) where {T, O <: Ring}
	make_αs(m.resolution.nside)
end
function make_αs(mm::MaskedMap{T,O}) where {T, O <: Ring}
	make_αs(mm.resolution.nside)[:,mm.inds]
end

# function make_αs(mm:MaskedMap{T,O}, coeff) where {T, O <: Ring}
# 	coeff .* healpy.alm2map_der1(healpy.synalm(
# 		read_Cℓs()), m.resolution.nside)[2:end, mm.inds]
	# end

###########################################################################
# Methods for plotting Masked Maps

plotmap(CMB::Map, c=:balance) = plot(CMB, c=c)
function GaussianLensing.plotmap(mm::MaskedMap; rot=0, center=nothing,
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
