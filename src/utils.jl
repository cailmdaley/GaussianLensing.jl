
function read_Cℓs(dir="/home/cailmdaley/spt/datafiles",
				  filename="CAMB_fiducial_cosmo_scalCls.dat")
	table = readdlm(joinpath(dir, filename))
	ℓmin = Int(table[:, 1][1])
	Dℓs = [zeros(ℓmin); table[:, 2]]
	ℓs   = 0:length(Dℓs)-1; 
	Cℓs = @. 2π/(ℓs * (ℓs + 1)) * Dℓs; Cℓs[1] = 0
	return Cℓs
end

#https://camb.readthedocs.io/en/latest/CAMBdemo.html
